# Validation of the false-positive rate on measured, signal-free spectra.
#
# The self-calibration of run_pipeline.jl uses simulated spectra and can only test the method under its
# own assumptions. The decisive test for an instrument is a set of real spectra that contain no signal
# (e.g. spectra taken with the signal path blocked, off-resonance configurations, or frequency ranges
# already excluded): analysed with exactly the settings of the search, their p values must be flat,
# with f(p <= 0.05) = 0.05 within its binomial error. This script runs the analysis chain on every CSV
# file of a directory and reports that.
#
# usage: julia --project=. -t auto scripts/validate_nulls.jl dir=path/to/spectra [options]
#   options as in run_pipeline.jl: n m kind tau K W estimator maxit tol pk_order
#   pattern=.csv   file suffix to select (default .csv)
#   out=...        per-spectrum results (CSV); default <output_root>/data/null_validation.csv
#   output_root=results

include(joinpath(@__DIR__, "run_pipeline.jl"))

const VALIDATE_USAGE = """
usage: julia --project=. -t auto scripts/validate_nulls.jl dir=path/to/spectra [pattern=.csv] [out=file.csv]
                    [n=4 m=64 kind=SG tau=1e-7 K=6 W=96 estimator=window maxit=60 tol=1e-10 pk_order=2] [output_root=results]
Runs the analysis chain on every signal-free spectrum in dir and tests that the p values are flat.
Use exactly the settings of the search. Files must be independent measurements.
"""

function validate_main(args)
    if any(a -> a in ("-h", "--help", "help"), args)
        print(VALIDATE_USAGE)
        return
    end
    kw = Dict{String,String}()
    for arg in args
        occursin('=', arg) || throw(ArgumentError("expected key=value, got '$arg'"))
        k, v = split(arg, '='; limit = 2)
        kw[k] = v
    end
    allowed = ("dir", "pattern", "out", "output_root", "n", "m", "kind", "tau", "K", "W", "estimator", "maxit", "tol", "pk_order")
    for key in keys(kw)
        key in allowed || throw(ArgumentError("unknown option '$key'"))
    end
    haskey(kw, "dir") || throw(ArgumentError("dir= is required\n" * VALIDATE_USAGE))
    dir = kw["dir"]; isdir(dir) || throw(ArgumentError("directory not found: $dir"))
    get(k, d) = haskey(kw, k) ? kw[k] : d
    paths = prepare_output_layout(get("output_root", Base.get(ENV, "MADMAX_OUT", "results")))
    open_log(joinpath(paths.logs, "validate_nulls.log"))
    try
        configure_threads()
        pattern = get("pattern", ".csv")
        files = sort([joinpath(dir, f) for f in readdir(dir) if endswith(f, pattern)])
        isempty(files) && throw(ArgumentError("no files ending in '$pattern' in $dir"))
        log_message("validate_nulls: $(length(files)) files in $dir")
        K = parse(Int, get("K", "6"))
        base = (; n = parse(Int, get("n", "4")), m = parse(Int, get("m", "64")), kind = Symbol(get("kind", "SG")),
                tau = parse(Float64, get("tau", "1e-7")), W = parse(Int, get("W", "96")),
                estimator = Symbol(get("estimator", "window")), boundary_removal = K > 0, nboundary = max(K, 1),
                pk_order = parse(Int, get("pk_order", "2")), maxit = parse(Int, get("maxit", "60")),
                tol = parse(Float64, get("tol", "1e-10")))
        # read all spectra; a file that cannot be read is reported and skipped (it is not a test result)
        spectra = Vector{Any}(undef, length(files))
        for (i, f) in enumerate(files)
            spectra[i] = try
                read_spectrum(f)
            catch err
                err isa InterruptException && rethrow()
                log_message("validate_nulls: skipping $f: $(sprint(showerror, err))")
                nothing
            end
        end
        # one filter setup per frequency grid
        setups = Dict{Tuple{Int,Float64,Float64},Any}()
        for sp in spectra
            sp === nothing && continue
            g = sp[2]; key = (g.nbins, g.df_hz, g.fref_hz)
            haskey(setups, key) && continue
            cfg = PipelineConfig(; base..., grid = g)
            setups[key] = filter_setup(cfg)
        end
        log_message(@sprintf("validate_nulls: %d frequency grid(s); SG n = %d, m = %d, tau = %.1e, K = %d",
                             length(setups), base.n, base.m, base.tau, K))
        status = fill("unread", length(files)); T = fill(NaN, length(files)); p = fill(NaN, length(files))
        s0 = fill(NaN, length(files)); s1 = fill(NaN, length(files)); f1 = fill(NaN, length(files)); f2 = fill(NaN, length(files))
        Threads.@threads for i in eachindex(files)
            sp = spectra[i]; sp === nothing && continue
            x, g = sp
            S = setups[(g.nbins, g.df_hz, g.fref_hz)]
            try
                r = run_chain(S, x)
                T[i] = r.T; p[i] = r.p_pk; s0[i] = r.noise.sigma0; s1[i] = r.noise.sigma1
                f1[i], f2[i] = run_frequencies(S.cfg, r.run)
                status[i] = "ok"
            catch err
                is_fit_failure(err) || rethrow()
                status[i] = "fit failed: " * replace(sprint(showerror, err), ',' => ';', '\n' => ' ')
            end
        end
        nread = count(!=("unread"), status); nok = count(==("ok"), status); nfail = nread - nok
        nread == 0 && error("no spectrum could be read")
        pv = p[status .== "ok"]
        # failed fits count as positives (conservative)
        f05 = (count(<=(0.05), pv) + nfail) / nread
        f01 = (count(<=(0.01), pv) + nfail) / nread
        err05 = sqrt(0.05 * 0.95 / nread)
        ks = ks_uniform(pv)
        hist = MR.p_histogram(pv)
        log_message(@sprintf("validate_nulls: %d spectra analysed, %d noise fits failed (counted as positives), %d files unreadable",
                             nok, nfail, length(files) - nread))
        log_message(@sprintf("validate_nulls: f(p <= 0.05) = %.4f ± %.4f (expected 0.05), f(p <= 0.01) = %.4f", f05, err05, f01))
        log_message(@sprintf("validate_nulls: KS distance from flat = %.4f, KS p = %.3g; p histogram (20 bins) %s", ks.D, ks.p, join(hist, " ")))
        nread < 100 && log_message(@sprintf("validate_nulls: only %d spectra: the false-positive rate is known to ±%.3f only; use at least a few hundred", nread, 3err05))
        ok = abs(f05 - 0.05) <= 3err05 && ks.p > 0.01
        log_message(ok ? "validate_nulls: PASSED, the false-positive rate on these signal-free instrument spectra agrees with the nominal one" :
                         "validate_nulls: FAILED, the p values on signal-free instrument spectra are not flat; the analysis is not calibrated for this instrument")
        out = get("out", joinpath(paths.data, "null_validation.csv"))
        MR.write_csv(out, ["file", "status", "T", "p_PK", "sigma0", "sigma1", "run_start_hz", "run_end_hz"],
                     hcat(files, status, T, p, s0, s1, f1, f2))
        log_message("validate_nulls: per-spectrum results written to $out")
        return (; files, status, T, p, f05, f01, err05, ks_D = ks.D, ks_p = ks.p, hist, passed = ok)
    finally
        close_log()
    end
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    validate_main(ARGS)
end
