function study_pvalues(o::StudyOptions)
    tmin, t0h, t0x = thresholds(o)
    cfg = base_config(o; boundary_removal = true)
    F = SGFilter(o.n, o.m; kind = o.kind); Dspk = spike_fourier(F, cfg.grid.nbins)
    S0 = filter_setup(with_tau(cfg, 0.01); with_pk = false, with_boundary = false)
    sig = axion_signal(); bmed = median_background()
    out = Dict{Tuple{String,Float64},Any}()
    for model in o.models
        taus = [tmin, model == Homoscedastic ? t0h : t0x, 0.01, 0.03, 0.05, 0.1, 0.2]
        setups = [filter_setup(with_tau(cfg, tau); Dspk) for tau in taus]
        for lv in o.levels
            cid = case_id(model, lv)
            bg = fit_set(S0, "noise", data_generator(o, cfg, cid; model, level = lv), o.nspec)
            sg = fit_set(S0, "signal", data_generator(o, cfg, cid; model, level = lv, signal = sig), o.nspec)
            b0 = fit_set(S0, "b0", data_generator(o, cfg, case_id(model, lv; background_on = false);
                                                  model, level = lv, background_on = false), o.nspec)
            gfilt = filt_generator(bg.sigma, o.m, o.seed + 17 + Int(model))
            fl = o.filt_reestimate ? fit_set(S0, "filt", gfilt, o.nfilt) :
                 SpectrumSet("filt", gfilt, [bg.sigma[mod1(i, length(bg))] for i in 1:o.nfilt])
            raw = map((bg, sg)) do set
                T = zeros(length(set))
                Threads.@threads for i in 1:length(set)
                    sigma = Float64.(set.sigma[i])
                    T[i] = squares_T(linear_residual(S0.F, set.gen(i))[S0.rows] ./ (sigma[S0.rows] .* sqrt(S0.q)))
                end
                pvalue.(Ref(setups[1].pbctab), T)
            end
            _, bref, sigmaref, sref = simulate_spectrum(spectrum_rng(o, cid, 1), cfg; model, level = lv, bmed, signal = sig)
            pertau = Dict{Float64,Any}()
            for (k, S) in enumerate(setups)
                iid = null_iid(S, o.nnull; seed = o.seed + 11)
                Tf, Tfp = eval_set(S, fl)
                nf = NullSet(:filt, sort(Tf)); nfp = NullSet(:filt, sort(Tfp))
                h = Dict{String,Any}()
                for (name, set) in (("noise", bg), ("signal", sg), ("b0", b0))
                    T, Tp = eval_set(S, set)
                    pf = mc_pvalue.(T, Ref(nf)); pfr = mc_pvalue.(Tp, Ref(nfp))
                    pii = mc_pvalue.(T, Ref(iid)); pir = pvalue.(Ref(S.pktab), Tp)
                    h[name] = (; filt = p_histogram(pf), filt_rm = p_histogram(pfr), iid = p_histogram(pii), pk_rm = p_histogram(pir),
                               f05 = (filt = fraction_below(pf, 0.05), filt_rm = fraction_below(pfr, 0.05),
                                      iid = fraction_below(pii, 0.05), pk_rm = fraction_below(pir, 0.05)))
                end
                lw = whitened_line(S, sref, sigmaref[o.m+1:end-o.m])
                c = argmax(abs.(lw))
                pertau[taus[k]] = (; hist = h, line = lw[max(1, c - 60):min(length(lw), c + 60)])
                log_message(@sprintf("pvalues %s %.2f%% \u03c4 = %.2e: f(p<0.05) noise filt %.3f iid %.3f | signal filt %.3f iid %.3f",
                                     model_tag(model), 100lv, taus[k], h["noise"].f05.filt, h["noise"].f05.iid,
                                     h["signal"].f05.filt, h["signal"].f05.iid))
            end
            out[(model_tag(model), lv)] = (; taus, pertau, raw_noise = p_histogram(raw[1]), raw_signal = p_histogram(raw[2]),
                                           raw_f95 = fraction_below(1 .- raw[1], 0.05))
        end
    end
    res = (; cases = out, tau_min = tmin, tau0 = (homo = t0h, hetero = t0x), levels = o.levels,
           models = model_tag.(o.models), nspec = o.nspec, nnull = o.nnull, nfilt = o.nfilt, A = sig[2])
    save_result(result_path(o, "pvalues"), res)
    return res
end
