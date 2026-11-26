function study_sigmafit(o::StudyOptions)
    cfg = base_config(o); F = SGFilter(o.n, o.m; kind = o.kind)
    wsg = WSGCache(F, cfg.grid.nbins); q = noise_pass_fraction(F); m = o.m
    bmed = median_background(); out = Dict{Tuple{String,Float64},Any}()
    for model in o.models, lv in o.levels
        x, b, sigma, _ = simulate_spectrum(spectrum_rng(o, case_id(model, lv), 1), cfg; model, level = lv, bmed)
        fw = estimate_noise(F, wsg, x; W = 96, estimator = :window, tol = cfg.tol, maxit = cfg.maxit)
        fl = estimate_noise(F, wsg, x; W = 96, estimator = :lsq, tol = cfg.tol, maxit = cfg.maxit)
        wt = 1 ./ pad_edges(fw.sigma, m) .^ 2
        _, _, sw96, bw96 = window_fit(x, 96, wt)
        _, _, sw256, bw256 = window_fit(x, 256, wt)
        r = x[m+1:end-m] .- fw.bhat
        nb = div(length(r), 96)
        brms = [sqrt(mean(abs2, view(r, (k-1)*96+1:k*96)) / q) for k in 1:nb]
        bblk = [mean(view(fw.bhat, (k-1)*96+1:k*96)) for k in 1:nb]
        sigma0, sigma1 = noise_parameters(model, lv, bmed)
        out[(model_tag(model), lv)] = (; sw96, bw96, sw256, bw256, brms, bblk,
            window = (fw.sigma0, fw.sigma1), lsq = (fl.sigma0, fl.sigma1), truth = (sigma0, sigma1),
            iterations = (fw.iterations, fl.iterations))
        log_message(@sprintf("sigmafit %s %.2f%%: window \u03c3\u03020 = %.4g \u03c3\u03021 = %.4g, lsq \u03c3\u03020 = %.4g \u03c3\u03021 = %.4g, truth %.4g %.4g",
                             model_tag(model), 100lv, fw.sigma0, fw.sigma1, fl.sigma0, fl.sigma1, sigma0, sigma1))
    end
    res = (; cases = out, bmed, levels = o.levels, models = model_tag.(o.models))
    save_result(result_path(o, "sigmafit"), res)
    return res
end

function study_iteration(o::StudyOptions; nsteps::Integer = 12)
    cfg = base_config(o); F = SGFilter(o.n, o.m; kind = o.kind)
    wsg = WSGCache(F, cfg.grid.nbins); bmed = median_background()
    ns = min(o.nspec, 100); out = Dict{Tuple{Symbol,Float64},Any}()
    for lv in o.levels, est in (:window, :lsq)
        d0 = fill(NaN, ns, nsteps); d1 = fill(NaN, ns, nsteps); conv = zeros(Int, ns)
        Threads.@threads for s in 1:ns
            x, = simulate_spectrum(spectrum_rng(o, case_id(Heteroscedastic, lv), s), cfg;
                                   model = Heteroscedastic, level = lv, bmed)
            nf = estimate_noise(F, wsg, x; W = 96, estimator = est, tol = 0.0, maxit = nsteps)
            k = length(nf.sigma0_hist) - 1
            d0[s, 1:k] = abs.(diff(nf.sigma0_hist)); d1[s, 1:k] = abs.(diff(nf.sigma1_hist))
            conv[s] = estimate_noise(F, wsg, x; W = 96, estimator = est, tol = cfg.tol, maxit = cfg.maxit).iterations
        end
        med(v) = (u = filter(isfinite, v); isempty(u) ? NaN : median(u))
        out[(est, lv)] = (; dsigma0 = [med(d0[:, k]) for k in 1:nsteps], dsigma1 = [med(d1[:, k]) for k in 1:nsteps],
                          iterations = conv)
        log_message(@sprintf("iteration %s %.2f%%: median iterations to tol = %.1f", est, 100lv, median(conv)))
    end
    res = (; cases = out, levels = o.levels, nsteps, nspec = ns)
    save_result(result_path(o, "iteration"), res)
    return res
end
