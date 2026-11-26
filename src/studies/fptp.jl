function study_fptp(o::StudyOptions)
    taus = tau_grid(o); cfg = base_config(o; boundary_removal = true); lv = o.level
    F = SGFilter(o.n, o.m; kind = o.kind); Dspk = spike_fourier(F, cfg.grid.nbins)
    S0 = filter_setup(with_tau(cfg, taus[1]); with_pk = false, with_boundary = false)
    sig = axion_signal()
    out = Dict{String,Any}()
    for model in o.models
        cid = case_id(model, lv)
        bg = fit_set(S0, "noise", data_generator(o, cfg, cid; model, level = lv), o.nspec)
        sg = fit_set(S0, "signal", data_generator(o, cfg, cid; model, level = lv, signal = sig), o.nspec)
        gfilt = filt_generator(bg.sigma, o.m, o.seed + 13 + Int(model))
        fl = o.filt_reestimate ? fit_set(S0, "filt", gfilt, o.nfilt) :
             SpectrumSet("filt", gfilt, [bg.sigma[mod1(i, length(bg))] for i in 1:o.nfilt])
        keys = ("iid_mc", "pk", "pk_rm", "filt_mc", "filt_mc_rm")
        FP = Dict(k => zeros(length(taus)) for k in keys); TP = Dict(k => zeros(length(taus)) for k in keys)
        for (k, tau) in enumerate(taus)
            S = filter_setup(with_tau(cfg, tau); Dspk)
            cutoff = decision_threshold(S.pktab, 0.05)
            iid = null_iid(S, o.nnull; seed = o.seed + 11)
            Tf, Tfp = eval_set(S, fl)
            nf = NullSet(:filt, sort(Tf)); nfp = NullSet(:filt, sort(Tfp))
            for (set, D) in ((bg, FP), (sg, TP))
                T, Tp = eval_set(S, set)
                D["iid_mc"][k] = mean(mc_pvalue.(T, Ref(iid)) .< 0.05)
                D["pk"][k] = mean(T .> cutoff)
                D["pk_rm"][k] = mean(Tp .> cutoff)
                D["filt_mc"][k] = mean(mc_pvalue.(T, Ref(nf)) .< 0.05)
                D["filt_mc_rm"][k] = mean(mc_pvalue.(Tp, Ref(nfp)) .< 0.05)
            end
            log_message(@sprintf("fptp %s \u03c4 = %.2e: FP iid %.3f pk %.3f pk_rm %.3f filt %.3f filt_rm %.3f | TP pk_rm %.3f",
                                 model_tag(model), tau, FP["iid_mc"][k], FP["pk"][k], FP["pk_rm"][k],
                                 FP["filt_mc"][k], FP["filt_mc_rm"][k], TP["pk_rm"][k]))
        end
        out[model_tag(model)] = (; FP, TP)
    end
    tmin = minimum_threshold(o)
    res = (; taus, cases = out, level = lv, A = sig[2], w = sig[3], tau_min = tmin, nspec = o.nspec,
           nnull = o.nnull, nfilt = o.nfilt, models = model_tag.(o.models))
    save_result(result_path(o, "fptp"), res)
    return res
end
