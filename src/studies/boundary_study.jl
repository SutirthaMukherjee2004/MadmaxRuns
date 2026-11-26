function study_boundary(o::StudyOptions; tau::Real = 1e-6)
    cfg = base_config(o; tau, boundary_removal = true)
    S = filter_setup(cfg); bmed = median_background(); lv = o.level
    x, = simulate_spectrum(spectrum_rng(o, case_id(Homoscedastic, lv), 1), cfg; model = Homoscedastic, level = lv, bmed)
    res = run_chain(S, x)
    B = res.boundary
    comp = res.z .- res.zp
    log_message(@sprintf("boundary \u03c4 = %.0e: s = %s; T = %.2f \u2192 %.2f", tau,
                         join([@sprintf("%.4g", v) for v in B.s[1:6]], ", "), squares_T(res.z), res.T))
    check = Dict{Int,Float64}()
    for K in (2, 6)
        SK = filter_setup(with(cfg; nboundary = K))
        p = zeros(o.nspec)
        Threads.@threads for i in 1:o.nspec
            xh, = simulate_spectrum(spectrum_rng(o, case_id(Heteroscedastic, lv), i), cfg; model = Heteroscedastic, level = lv, bmed)
            p[i] = run_chain(SK, xh).p_pk
        end
        check[K] = fraction_below(p, 0.05)
        log_message(@sprintf("boundary hetero K = %d: f(p<0.05) = %.4f", K, check[K]))
    end
    out = (; tau, z = res.z, component = comp, zp = res.zp, T_before = squares_T(res.z), T_after = res.T,
           s = B.s, uncaptured = B.uncaptured, patterns = B.V, fp_K = check, rows = collect(S.rows))
    save_result(result_path(o, "boundary"), out)
    return out
end

function study_tau0(o::StudyOptions; ns = 2:2:16, ms = collect(16:8:200), alpha::Real = 0.05)
    bmed = median_background(); lv = o.level
    xs = Dict(model => simulate_spectrum(spectrum_rng(o, case_id(model, lv), 1), base_config(o); model, level = lv, bmed)[1]
              for model in (Homoscedastic, Heteroscedastic))
    function fitted_sigma(F, model)
        nf = estimate_noise(F, WSGCache(F, length(xs[model])), xs[model])
        return pad_edges(nf.sigma, F.m)
    end
    tau0 = Dict{Tuple{Int,Int,String},Float64}()
    map0 = Dict(t => fill(NaN, length(ns), length(ms)) for t in ("homo", "hetero"))
    for (a, n) in enumerate(ns), (j, m) in enumerate(ms)
        2m + 1 > n + 1 || continue
        F = SGFilter(n, m; kind = o.kind); N = PAPER_GRID.nbins - 2m
        D = spike_fourier(F, PAPER_GRID.nbins); G = transfer_function(F, N)
        for model in (Homoscedastic, Heteroscedastic)
            t = tau0_exact(F, D, G, fitted_sigma(F, model); alpha)
            map0[model_tag(model)][a, j] = t; tau0[(n, m, model_tag(model))] = t
        end
    end
    F = SGFilter(o.n, o.m; kind = o.kind); N = PAPER_GRID.nbins - 2o.m
    D = spike_fourier(F, PAPER_GRID.nbins); G = transfer_function(F, N)
    tgrid = exp10.(range(-7, -1.5, length = 90))
    cross = Dict{String,Any}()
    for model in (Homoscedastic, Heteroscedastic)
        sigmaf = fitted_sigma(F, model)
        sig = zeros(length(tgrid)); zalpha = zeros(length(tgrid))
        for k in eachindex(tgrid)
            sig[k], zalpha[k], _, _ = boundary_significance(F, D, G, sigmaf, tgrid[k]; alpha)
        end
        t0 = tau0_exact(F, D, G, sigmaf; alpha)
        tau0[(o.n, o.m, model_tag(model))] = t0
        r2 = boundary_r2(sigmaf[o.m+1], sigmaf[end-o.m])
        closed = haskey(TAU0_CONSTANTS, o.n) ? (TAU0_CONSTANTS[o.n][1] * sqrt(r2)) .* tgrid .^ (-TAU0_CONSTANTS[o.n][2]) : fill(NaN, length(tgrid))
        tc = haskey(TAU0_CONSTANTS, o.n) ? tau0_closed(o.n, o.m, N - 2o.m, r2; alpha) : NaN
        cross[model_tag(model)] = (; significance = sig, closed, z_alpha = zalpha, tau0 = t0, tau0_closed = tc, r2)
        log_message(@sprintf("tau0 %s (n = %d, m = %d): exact %.2e, closed form %.2e", model_tag(model), o.n, o.m, t0, tc))
    end
    tmin, = thresholds(o)
    res = (; ns = collect(ns), ms, map = map0, tau0, taus = tgrid, cross, n = o.n, m = o.m, tau_min = tmin, alpha = alpha)
    save_result(result_path(o, "tau0"), res)
    return res
end
