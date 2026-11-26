function leakage_fraction(S::FilterSetup, b::AbstractVector, sigma::AbstractVector)
    C = whitened_fourier(rfft(linear_residual(S.F, b)), S.G, S.keep)
    c = (irfft(C, S.N) ./ sigma)[S.rows]
    if S.cfg.boundary_removal
        B = boundary_patterns(S.H, pad_edges(sigma, S.cfg.m)[S.J], sigma[S.rows], S.kappa; K = S.cfg.nboundary)
        c = remove_boundary(c ./ sqrt(S.kappa), B) .* sqrt(S.kappa)
    end
    v = var(c)
    return v / (v + S.kappa)
end

function study_leakage(o::StudyOptions; ns = 2:8, ms = collect(8:8:296), nbkg::Integer = 20)
    lv = minimum(o.levels); bmed = median_background()
    sigma0, sigma1 = noise_parameters(Homoscedastic, lv, bmed)
    bk = nbkg == 1 ? [nominal_background()] : [nominal_background(), (background(spectrum_rng(o, 950, k)) for k in 2:nbkg)...]
    out = Dict{Symbol,Matrix{Float64}}()
    for kind in (:SG, :SGW)
        phi = fill(NaN, length(ns), length(ms))
        for (a, n) in enumerate(ns), (j, m) in enumerate(ms)
            2m + 1 > n + 1 || continue
            cfg = PipelineConfig(n = n, m = m, kind = kind, tau = 1e-7, boundary_removal = true)
            S = filter_setup(cfg; with_pk = false, with_boundary = true)
            phi[a, j] = mean(leakage_fraction(S, b, fill(sigma0, S.N) .+ sigma1 .* b[m+1:end-m]) for b in bk)
        end
        out[kind] = phi
        log_message("leakage $kind done")
    end
    res = (; ns = collect(ns), ms, phi = out, level = lv, nbkg)
    save_result(result_path(o, "leakage"), res)
    return res
end
