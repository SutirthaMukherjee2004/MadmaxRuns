using MadmaxRuns, DelimitedFiles, Random, Printf, Statistics, LinearAlgebra
const MR = MadmaxRuns
include(joinpath(@__DIR__, "spectrum_calibration.jl"))

function read_spectrum(path)
    isfile(path) || throw(ArgumentError("spectrum CSV does not exist: $path"))
    d = try
        readdlm(path, ',', Float64; comments = true)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("could not read $path as a headerless numeric CSV: $(sprint(showerror, err))"))
    end
    size(d, 2) in (1, 2) || throw(ArgumentError("CSV must contain power or frequency_hz,power"))
    all(v -> v isa Real && isfinite(v), d) || throw(ArgumentError("CSV must contain only finite numbers"))
    if size(d, 2) == 1
        size(d, 1) == MR.PAPER_GRID.nbins ||
            throw(ArgumentError("power-only CSV must have $(MR.PAPER_GRID.nbins) rows; use frequency_hz,power for other grids"))
        return vec(d), MR.PAPER_GRID
    end
    size(d, 1) >= 2 || throw(ArgumentError("CSV needs at least two frequency bins"))
    f = d[:, 1]; x = d[:, 2]
    spacings = diff(f)
    df = first(spacings)
    tol = max(1e-8 * df, 8eps(maximum(abs, f)))
    all(s -> isfinite(s) && s > 0 && abs(s - df) <= tol, spacings) ||
        throw(ArgumentError("frequencies must be strictly increasing and uniformly spaced in Hz"))
    return x, MR.SpectrumGrid(length(x), df, f[1])
end

function bool_option(value)
    lowercase(value) in ("true", "1", "yes") && return true
    lowercase(value) in ("false", "0", "no") && return false
    throw(ArgumentError("expected true or false, got '$value'"))
end

function process_spectrum(kw, checkpoint_dir, data_dir)
    get(k, d) = haskey(kw, k) ? kw[k] : d
    n = parse(Int, get("n", "4")); m = parse(Int, get("m", "64"))
    K = parse(Int, get("K", "6")); seed = parse(Int, get("seed", "2026"))
    nnull = parse(Int, get("nnull", "0")); nfilt = parse(Int, get("nfilt", "0"))
    nnull >= 0 && nfilt >= 0 || throw(ArgumentError("null replica counts must be nonnegative"))
    K >= 0 || throw(ArgumentError("K must be nonnegative"))
    filt_reestimate = bool_option(get("filt_reestimate", "false"))
    resume = bool_option(get("resume", "true"))
    calibrate = bool_option(get("calibrate", "true"))
    bkg_model = Symbol(lowercase(get("bkg_model", "paper")))
    bkg_model in (:paper, :data) || throw(ArgumentError("bkg_model must be paper or data"))
    calib_signal = parse(Float64, get("calib_signal", "3"))
    isfinite(calib_signal) && calib_signal >= 0 || throw(ArgumentError("calib_signal must be a nonnegative peak in units of sigma"))
    nbkg = !calibrate ? 0 : haskey(kw, "nbkg") ? parse(Int, kw["nbkg"]) : ask_background_count(5000)
    calibrate && nbkg <= 0 && throw(ArgumentError("nbkg must be positive"))
    check = bool_option(get("check", "false"))
    tauscan = bool_option(get("tauscan", "false"))
    taus = haskey(kw, "taus") ? sort(unique(parse.(Float64, split(kw["taus"], ',')))) : tauscan ? copy(DEFAULT_SCAN_TAUS) : Float64[]
    all(t -> isfinite(t) && 0 < t < 1, taus) || throw(ArgumentError("taus must lie in (0, 1)"))
    check && !calibrate && throw(ArgumentError("check=true compares the spectrum with the noise-only mocks and needs calibrate=true"))
    configure_threads()
    x, grid = haskey(kw, "input") ? read_spectrum(kw["input"]) : (nothing, MR.PAPER_GRID)
    cfg = PipelineConfig(n = n, m = m, kind = Symbol(get("kind", "SG")), tau = parse(Float64, get("tau", "1e-7")),
                         W = parse(Int, get("W", "96")), estimator = Symbol(get("estimator", "window")),
                         boundary_removal = K > 0, nboundary = max(K, 1), pk_order = parse(Int, get("pk_order", "2")),
                         maxit = parse(Int, get("maxit", "60")), tol = parse(Float64, get("tol", "1e-10")), grid = grid)
    grid.nbins > 4m || throw(ArgumentError("the spectrum must have more than 4m bins"))
    cfg.estimator in (:window, :lsq) || throw(ArgumentError("estimator must be window or lsq"))
    cfg.estimator != :window || 4 <= cfg.W <= div(grid.nbins, 2) ||
        throw(ArgumentError("window fit needs W >= 4 and at least two full windows"))
    cfg.maxit > 0 && isfinite(cfg.tol) && cfg.tol >= 0 ||
        throw(ArgumentError("maxit must be positive and tolerance finite and nonnegative"))
    if x === nothing
        model = parse_noise_model(get("model", "homo")); level = parse(Float64, get("level", "0.05"))
        peak = parse(Float64, get("signal", "0"))
        isfinite(level) && level > 0 && isfinite(peak) && peak >= 0 ||
            throw(ArgumentError("mock noise level must be positive and signal peak nonnegative"))
        sig = peak > 0 ? (MR.AxionLine(), line_area_for_peak(MR.AxionLine(), 3.4, peak), 3.4, MR.axion_offset(grid)) : nothing
        x, = simulate_spectrum(MR.stream_rng(seed, 1), cfg; model, level, signal = sig)
        log_message(@sprintf("step 1: simulated spectrum (%s, \u03c3/b_med = %.2f %%, axion peak %.1f \u03c3)", model, 100level, peak))
    end
    all(isfinite, x) || throw(ArgumentError("spectrum contains non-finite values"))
    calib = (; nbkg, bkg_model, peak = calib_signal, check, taus)
    run = MR.CheckpointRun(checkpoint_dir, (; cfg, x, seed, nnull, nfilt, filt_reestimate, calib); resume)
    kind = haskey(kw, "input") ? "measured" : "mock"
    out = get("out", joinpath(data_dir, "$(kind)_$(run.signature[1:12]).jls"))
    log_message("Run checkpoint directory: $(run.dir)")
    return MR.with_checkpoint_lock(joinpath(run.dir, "pipeline.lock")) do
        score_spectrum(run, cfg, x; seed, nnull, nfilt, filt_reestimate, out, calib)
    end
end

function checkpointed_null(f, run, kind, n)
    statistics = Float64[]
    sizehint!(statistics, n)
    for first_replica in 1:128:n
        count = min(128, n - first_replica + 1)
        block = MR.run_stage(run, Symbol(kind, "_null_", first_replica)) do
            f(count, first_replica)
        end
        append!(statistics, block.T)
    end
    return NullSet(kind, sort!(statistics))
end

function score_spectrum(run, cfg, x; seed, nnull, nfilt, filt_reestimate, out, calib)
    input = MR.run_stage(run, :input) do
        (; cfg, x)
    end
    x = input.x
    S = MR.run_stage(run, :filter_setup) do
        filter_setup(cfg; with_pk = false)
    end
    log_message(@sprintf("setup: %s filter, n = %d, m = %d, N = %d, N_g = %d, q = %.5f, \u03c4 = %.2e, %d rfft modes removed, \u03ba = %.5f",
                         cfg.kind, cfg.n, cfg.m, S.N, S.Ng, S.q, cfg.tau, MR.n_removed_rfft(S.keep), S.kappa))
    r = MR._run_chain(S, x, run)
    nf = r.noise
    log_message(@sprintf("steps 2\u20133: \u03c3\u03020 = %.6g, \u03c3\u03021 = %.6g after %d iterations (converged: %s, estimator %s)",
                         nf.sigma0, nf.sigma1, nf.iterations, nf.converged, nf.estimator))
    log_message("          \u03c3\u03020 history: " * join([@sprintf("%.8g", v) for v in nf.sigma0_hist], " "))
    log_message("          \u03c3\u03021 history: " * join([@sprintf("%.8g", v) for v in nf.sigma1_hist], " "))
    if r.boundary !== nothing
        log_message("step 8: boundary singular values " * join([@sprintf("%.4g", v) for v in r.boundary.s[1:min(6, end)]], ", "))
    end
    run_freq = run_frequencies(cfg, r.run)
    log_message(@sprintf("step 9: T = %.4f (RunStatistics.t_obs), run %d\u2013%d of the tested bins = %.6f\u2013%.6f MHz (%d bins)",
                         r.T, r.run..., run_freq[1] / 1e6, run_freq[2] / 1e6, r.run[2] - r.run[1] + 1))
    log10p = pk_log10pvalue(S.pk, r.T)
    beyond = r.T > S.pk.tail_t
    log_message(@sprintf("step 10: p = %.5g (log10 p = %.3f) from P_K (order %d)%s; raw residual T = %.3f, p = %.4g from P_BC",
                         r.p_pk, log10p, cfg.pk_order,
                         beyond ? @sprintf(" [beyond the validated range p >= %.0e: conservative extrapolation, significance not calibrated at this level]", MR.PK_TAIL_SWITCH) : "",
                         r.T_raw, r.p_raw))
    p_iid = NaN; p_filt = NaN
    if nnull > 0
        iid = checkpointed_null(run, :iid, nnull) do count, first_replica
            null_iid(S, count; seed = seed + 11, first_replica)
        end
        p_iid = mc_pvalue(r.T, iid)
        log_message(@sprintf("H0^iid Monte Carlo (%d replicas): p = %.5g", nnull, p_iid))
    end
    if nfilt > 0
        filt = checkpointed_null(run, :filt, nfilt) do count, first_replica
            null_filt(S, [nf.sigma], count; seed = seed + 13, reestimate = filt_reestimate, first_replica)
        end
        p_filt = mc_pvalue(r.T, filt)
        log_message(@sprintf("H0^filt Monte Carlo (%d replicas): p = %.5g", nfilt, p_filt))
    end
    scan = isempty(calib.taus) ? nothing : scan_tau(S, cfg, x, nf, calib.taus)
    calibration = calib.nbkg > 0 ?
        calibrate_spectrum(run, S, cfg, nf, r; nbkg = calib.nbkg, bkg_model = calib.bkg_model, peak = calib.peak, seed,
                           check = calib.check, scan) :
        nothing
    if scan !== nothing && calibration === nothing
        log_message("tau scan: no calibration (calibrate=false): the smallest p over the scan is not corrected for the number of thresholds tried")
    end
    summary = result_summary(cfg, S, r, log10p, run_freq, calibration)
    if !isempty(out)
        log_message("Writing final result to $out")
        save_result(out, (; cfg, result = r, log10p, run_frequency_hz = run_freq, p_iid, p_filt, calibration, summary,
                          tauscan = scan === nothing ? nothing : scan.rows,
                          checkpoint_signature = run.signature))
    end
end

# Threshold scan for one spectrum and one SG filter (tauscan=true or taus=...). The noise fit does not
# depend on tau (steps 2-3), so it is done once and reused; steps 4-10 are repeated at every tau with
# that tau's Fourier mask, boundary patterns and P_K. With the boundary term removed, every small tau is
# valid (paper Sec. 3.7, Fig. 11) and the scan shows how stable a result is. Taking the smallest p
# over the scan is a choice among several correlated tests; the calibration therefore repeats the same
# scan on every noise-only realisation and reports p_scan, the trials-corrected probability.
const DEFAULT_SCAN_TAUS = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.05, 0.2]

function scan_tau(S, cfg, x, nf, taus)
    log_message(@sprintf("tau scan: %d thresholds %s, SG n = %d, m = %d, same noise fit", length(taus),
                         join([@sprintf("%.0e", t) for t in taus], ","), cfg.n, cfg.m))
    Dspk = cfg.boundary_removal ? S.Dspk : nothing
    setups = [filter_setup(MR.with_tau(cfg, t); with_boundary = cfg.boundary_removal, Dspk) for t in taus]
    rows = NamedTuple[]
    for (t, St) in zip(taus, setups)
        rr = run_chain(St, x; noise = nf)
        fr = run_frequencies(cfg, rr.run)
        lp = pk_log10pvalue(St.pk, rr.T)
        push!(rows, (; tau = t, removed = MR.n_removed_rfft(St.keep), kappa = St.kappa, T = rr.T, p = rr.p_pk, log10p = lp,
                       beyond = rr.T > St.pk.tail_t, run_start_hz = fr[1], run_end_hz = fr[2]))
        log_message(@sprintf("tau scan:   tau = %.1e  removed %4d  kappa %.4f  T = %8.3f  p = %.4g%s  run %.6f-%.6f MHz",
                             t, rows[end].removed, St.kappa, rr.T, rr.p_pk, rows[end].beyond ? " (beyond validated range)" : "",
                             fr[1] / 1e6, fr[2] / 1e6))
    end
    k = argmin([r.log10p for r in rows])
    log_message(@sprintf("tau scan: smallest p = %.4g at tau = %.1e (not trials-corrected; see p_scan from the calibration)", rows[k].p, rows[k].tau))
    return (; taus, setups, rows, min_log10p = rows[k].log10p, best_tau = rows[k].tau)
end

# First and last frequency (Hz) of the strongest run; tested bin j is row j + 2m of the input.
run_frequencies(cfg, run) = run[1] > 0 ?
    (cfg.grid.fref_hz + (run[1] + 2cfg.m - 1) * cfg.grid.df_hz, cfg.grid.fref_hz + (run[2] + 2cfg.m - 1) * cfg.grid.df_hz) : (NaN, NaN)

# Final summary: what was measured, which checks passed, and down to which p the significance is
# calibrated (P_K: validated to 1e-4 in App. G; p_bkg: resolution 1/(nbkg + 1)).
function result_summary(cfg, S, r, log10p, run_freq, cal)
    checked = cal !== nothing && cal.check !== nothing
    status = cal === nothing ? "NOT CALIBRATED (calibrate=false)" :
             !cal.calibrated ? "NOT TRUSTWORTHY: calibration failed" :
             checked && !cal.check.passed ? "NOT TRUSTWORTHY: the data violate the assumptions (check failed)" :
             !checked ? "CALIBRATED, assumptions not checked (run with check=true)" : "CALIBRATED AND CHECKED"
    floor_pk = MR.PK_TAIL_SWITCH
    floor_mc = cal === nothing ? NaN : 1 / (cal.nbkg + 1)
    log_message("RESULT: " * status)
    log_message(@sprintf("RESULT: T = %.4f, strongest run %.6f\u2013%.6f MHz, p(P_K) = %.4g (log10 %.2f)%s",
                         r.T, run_freq[1] / 1e6, run_freq[2] / 1e6, r.p_pk, log10p,
                         r.T > S.pk.tail_t ? ", below the validated range" : ""))
    cal === nothing || log_message(@sprintf("RESULT: p_bkg %s %.4g from %d background realisations",
                                            cal.p_bkg_limit ? "<" : "=", cal.p_bkg, cal.nbkg))
    cal === nothing || cal.scan === nothing ||
        log_message(@sprintf("RESULT: tau scan: smallest p %.4g at tau = %.1e; trials-corrected p_scan %s %.4g",
                             cal.scan.min_p, cal.scan.best_tau, cal.scan.p_scan_limit ? "<" : "=", cal.scan.p_scan))
    log_message(@sprintf("RESULT: significance is calibrated down to p = %.0e (P_K)%s; smaller p values are not calibrated by this run",
                         floor_pk, cal === nothing ? "" : @sprintf(" and p = %.1e (Monte Carlo)", floor_mc)))
    log_message("RESULT: before claiming a detection, validate the false-positive rate on signal-free instrument spectra (scripts/validate_nulls.jl) and correct for the number of spectra searched")
    return (; status, T = r.T, p_pk = r.p_pk, log10p, run_frequency_hz = run_freq,
            p_bkg = cal === nothing ? NaN : cal.p_bkg, p_bkg_limit = cal === nothing ? false : cal.p_bkg_limit,
            calibrated_to_pk = floor_pk, calibrated_to_mc = floor_mc)
end

const USAGE = """
usage: julia --project=. -t auto scripts/run_pipeline.jl [key=value ...]

Real spectrum:   input=spectrum.csv   (CSV: frequency_hz,power; one row per bin)
Mock spectrum:   model=homo|hetero level=0.05 signal=0   (no input=)

Analysis:     n=4 m=64 kind=SG|SGW tau=1e-7 K=6 W=96 estimator=window|lsq maxit=60 tol=1e-10 pk_order=2 seed=2026
Calibration:  calibrate=true nbkg=5000 (asked in a terminal if omitted) bkg_model=paper|data calib_signal=3
Check:        check=false   (true: test the measured residual against the noise-only mocks; recommended for real data)
Tau scan:     tauscan=false (true: repeat the analysis at tau = 1e-7..0.2 with the same SG filter) or taus=1e-7,1e-5,...
MC nulls:     nnull=0 nfilt=0 filt_reestimate=false
Output:       output_root=results out=... log=... checkpoint_dir=... resume=true

Example:      julia --project=. -t auto scripts/run_pipeline.jl input=spec.csv nbkg=5000 check=true
Full documentation: README.md
"""

function main(args)
    if any(a -> a in ("-h", "--help", "help"), args)
        print(USAGE)
        return
    end
    kw = Dict{String,String}()
    for arg in args
        occursin('=', arg) || throw(ArgumentError("expected key=value, got '$arg'"))
        k, v = split(arg, '='; limit = 2)
        haskey(kw, k) && throw(ArgumentError("duplicate option '$k'"))
        kw[k] = v
    end
    paths = prepare_output_layout(get(kw, "output_root", get(ENV, "MADMAX_OUT", "results")))
    checkpoint_dir = get(kw, "checkpoint_dir", joinpath(paths.data, "checkpoints"))
    logfile = get(kw, "log", joinpath(paths.logs, "pipeline.log"))
    open_log(logfile)
    try
        log_message("Starting single-spectrum run in $(paths.root)")
        allowed = ("input", "out", "n", "m", "K", "seed", "nnull", "nfilt", "kind", "tau", "W",
                   "estimator", "pk_order", "maxit", "tol", "model", "level", "signal",
                   "checkpoint_dir", "log", "resume", "filt_reestimate", "output_root",
                   "calibrate", "nbkg", "bkg_model", "calib_signal", "check", "tauscan", "taus")
        for key in keys(kw)
            key in allowed || throw(ArgumentError("unknown pipeline option '$key'"))
        end
        process_spectrum(kw, checkpoint_dir, paths.data)
    catch err
        log_message("Pipeline stopped: $(sprint(showerror, err))")
        out = get(kw, "out", "")
        isempty(out) || !isfile(out) || log_message("An existing output file at $out may belong to an earlier run")
        rethrow()
    finally
        close_log()
    end
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Base.exit_on_sigint(false)
    main(ARGS)
end
