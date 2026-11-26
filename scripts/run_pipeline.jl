using MadmaxRuns, DelimitedFiles, Random, Printf
const MR = MadmaxRuns

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
    run = MR.CheckpointRun(checkpoint_dir, (; cfg, x, seed, nnull, nfilt, filt_reestimate); resume)
    kind = haskey(kw, "input") ? "measured" : "mock"
    out = get("out", joinpath(data_dir, "$(kind)_$(run.signature[1:12]).jls"))
    log_message("Run checkpoint directory: $(run.dir)")
    return MR.with_checkpoint_lock(joinpath(run.dir, "pipeline.lock")) do
        score_spectrum(run, cfg, x; seed, nnull, nfilt, filt_reestimate, out)
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

function score_spectrum(run, cfg, x; seed, nnull, nfilt, filt_reestimate, out)
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
    log_message(@sprintf("step 9: T = %.4f (RunStatistics.t_obs), run %d\u2013%d of the tested bins", r.T, r.run...))
    log_message(@sprintf("step 10: p = %.5g from P_K (order %d); raw residual T = %.3f, p = %.4g from P_BC",
                         r.p_pk, cfg.pk_order, r.T_raw, r.p_raw))
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
    if !isempty(out)
        log_message("Writing final result to $out")
        save_result(out, (; cfg, result = r, p_iid, p_filt, checkpoint_signature = run.signature))
    end
end

function main(args)
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
                   "checkpoint_dir", "log", "resume", "filt_reestimate", "output_root")
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
