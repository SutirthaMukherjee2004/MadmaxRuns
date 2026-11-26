excess_kurtosis(v) = (mu = mean(v); mean((v .- mu) .^ 4) / mean((v .- mu) .^ 2)^2 - 3)

function study_examples(o::StudyOptions)
    cfg = base_config(o; boundary_removal = true)
    S = filter_setup(cfg); sig = axion_signal(); bmed = median_background(); lv = o.level
    log_message(@sprintf("examples: \u03c4 = %.1e, %d rfft modes removed, \u03ba = %.4f", o.tau, n_removed_rfft(S.keep), S.kappa))
    out = Dict{String,Any}()
    for model in o.models
        cid = case_id(model, lv)
        x, b, sigma, s = simulate_spectrum(spectrum_rng(o, cid, 1), cfg; model, level = lv, bmed, signal = sig)
        res = run_chain(S, x); nf = res.noise; m = o.m
        sigmabar = median(nf.sigma)
        pre = sigmabar .* res.r ./ nf.sigma
        wt = 1 ./ pad_edges(nf.sigma, m) .^ 2
        _, _, sw, bw = window_fit(x, cfg.W, wt)
        N = S.N; fr = collect(0:div(N, 2)) ./ N
        xin = x[m+1:end-m]
        one = (; x, b, sigma = sigma, line = s, bhat = nf.bhat, sigma_hat = nf.sigma, fit = (nf.sigma0, nf.sigma1),
               truth = noise_parameters(model, lv, bmed), f = fr,
               Px = abs2.(rfft(xin .- mean(xin))) ./ N, Pr = abs2.(res.R) ./ N, Pmodel = mean(abs2, nf.sigma) .* S.G .^ 2,
               r = res.r, pre = pre, kurt_r = excess_kurtosis(res.r), kurt_pre = excess_kurtosis(pre),
               sw, bw, before = abs2.(res.R) ./ (N * sigmabar^2), after = abs2.(res.C) ./ (N * sigmabar^2),
               removed = findall(.!S.keep) .- 1, T = res.T, p = res.p_pk, z = res.zp)
        p_raw0 = zeros(o.nspec); p0 = zeros(o.nspec); p1 = zeros(o.nspec)
        Threads.@threads for i in 1:o.nspec
            xs, = simulate_spectrum(spectrum_rng(o, cid, i), cfg; model, level = lv, bmed)
            r0 = run_chain(S, xs); p_raw0[i] = r0.p_raw; p0[i] = r0.p_pk
            xl, = simulate_spectrum(spectrum_rng(o, cid, i), cfg; model, level = lv, bmed, signal = sig)
            p1[i] = run_chain(S, xl).p_pk
        end
        out[model_tag(model)] = (; one..., p_raw = p_raw0, p_noise = p0, p_signal = p1,
                                 fp = fraction_below(p0, 0.05), tp = fraction_below(p1, 0.05),
                                 raw_above95 = count(>(0.95), p_raw0) / o.nspec)
        log_message(@sprintf("examples %s: \u03c3\u03020 = %.3f \u03c3\u03021 = %.4g; kurtosis r %.2f, pre-normalised %.2f; raw p>0.95 %.3f; FP %.3f TP %.3f",
                             model_tag(model), nf.sigma0, nf.sigma1, one.kurt_r, one.kurt_pre,
                             out[model_tag(model)].raw_above95, out[model_tag(model)].fp, out[model_tag(model)].tp))
    end
    res = (; cases = out, tau = o.tau, kappa = S.kappa, level = lv, nspec = o.nspec, A = sig[2])
    save_result(result_path(o, "examples"), res)
    return res
end
