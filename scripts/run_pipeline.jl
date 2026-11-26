using MadmaxRuns, DelimitedFiles, Random, Printf
const MR = MadmaxRuns

function read_spectrum(path)
    d = readdlm(path, ',', Float64; comments = true)
    size(d, 2) == 1 && return vec(d), MR.PAPER_GRID
    f = d[:, 1]; x = d[:, 2]
    return x, MR.SpectrumGrid(length(x), f[2] - f[1], f[1])
end

function main(args)
    kw = Dict{String,String}(String.(split(a, '='; limit = 2)) for a in args)
    get(k, d) = haskey(kw, k) ? kw[k] : d
    n = parse(Int, get("n", "4")); m = parse(Int, get("m", "64"))
    K = parse(Int, get("K", "6")); seed = parse(Int, get("seed", "2026"))
    configure_threads()
    x, grid = haskey(kw, "input") ? read_spectrum(kw["input"]) : (nothing, MR.PAPER_GRID)
    cfg = PipelineConfig(n = n, m = m, kind = Symbol(get("kind", "SG")), tau = parse(Float64, get("tau", "1e-7")),
                         W = parse(Int, get("W", "96")), estimator = Symbol(get("estimator", "window")),
                         boundary_removal = K > 0, nboundary = max(K, 1), pk_order = parse(Int, get("pk_order", "2")),
                         grid = grid)
    if x === nothing
        model = parse_noise_model(get("model", "homo")); level = parse(Float64, get("level", "0.05"))
        peak = parse(Float64, get("signal", "0"))
        sig = peak > 0 ? (MR.AxionLine(), line_area_for_peak(MR.AxionLine(), 3.4, peak), 3.4, MR.axion_offset(grid)) : nothing
        x, = simulate_spectrum(MR.stream_rng(seed, 1), cfg; model, level, signal = sig)
        log_message(@sprintf("step 1: simulated spectrum (%s, \u03c3/b_med = %.2f %%, axion peak %.1f \u03c3)", model, 100level, peak))
    end
    S = filter_setup(cfg)
    log_message(@sprintf("setup: %s filter, n = %d, m = %d, N = %d, N_g = %d, q = %.5f, \u03c4 = %.2e, %d rfft modes removed, \u03ba = %.5f",
                         cfg.kind, n, m, S.N, S.Ng, S.q, cfg.tau, MR.n_removed_rfft(S.keep), S.kappa))
    r = run_chain(S, x)
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
    nnull = parse(Int, get("nnull", "0")); nfilt = parse(Int, get("nfilt", "0"))
    p_iid = NaN; p_filt = NaN
    if nnull > 0
        p_iid = mc_pvalue(r.T, null_iid(S, nnull; seed = seed + 11))
        log_message(@sprintf("H0^iid Monte Carlo (%d replicas): p = %.5g", nnull, p_iid))
    end
    if nfilt > 0
        p_filt = mc_pvalue(r.T, null_filt(S, [nf.sigma], nfilt; seed = seed + 13))
        log_message(@sprintf("H0^filt Monte Carlo (%d replicas): p = %.5g", nfilt, p_filt))
    end
    out = get("out", "")
    isempty(out) || save_result(out, (; cfg, result = r, p_iid, p_filt))
end

main(ARGS)
