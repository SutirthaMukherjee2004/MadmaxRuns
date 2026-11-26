function study_background(o::StudyOptions)
    g = PAPER_GRID; f = collect(frequencies(g)); lv = 0.01
    rng = spectrum_rng(o, 900, 1)
    reals = [background(rng) for _ in 1:10]
    bnom = nominal_background(); bmed = median_background()
    sigma0, sigma1 = noise_parameters(Homoscedastic, lv, bmed)
    sigma = noise_sigma(sigma0, sigma1, bnom)
    x = bnom .+ sigma .* randn(rng, g.nbins)
    c = frequency_bin(g, axion_offset(g)) - 0.5
    A = line_area_for_peak(AxionLine(), 3.4, 3.0)
    lines = [sigma .* line_profile(AxionLine(), A, 3.4k, g.nbins, c) for k in (1, 4, 16)]
    F = SGFilter(o.n, o.m)
    r = linear_residual(F, x)
    h = o.m; env = [std(view(r, max(1, i - h):min(length(r), i + h))) for i in eachindex(r)]
    res = (; f, reals, nominal = bnom, x, sigma = sigma, lines, centre = c, widths = 3.4 .* [1, 4, 16], A,
           r, envelope = env, f_r = f[o.m+1:end-o.m], level = lv)
    save_result(result_path(o, "background"), res)
    return res
end

function study_spectra(o::StudyOptions)
    g = PAPER_GRID; bmed = median_background(); cfg = base_config(o)
    out = Dict{Tuple{String,Float64},Any}()
    for model in o.models, lv in o.levels
        x, b, sigma, _ = simulate_spectrum(spectrum_rng(o, case_id(model, lv), 1), cfg; model, level = lv, bmed)
        out[(model_tag(model), lv)] = (; x, b, sigma = sigma)
    end
    res = (; f = collect(frequencies(g)), spectra = out, levels = o.levels, models = model_tag.(o.models))
    save_result(result_path(o, "spectra"), res)
    return res
end
