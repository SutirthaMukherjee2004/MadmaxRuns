running_max(v, h) = [maximum(abs, view(v, max(1, i - h):min(length(v), i + h))) for i in eachindex(v)]
function running_mean(v, h)
    n = length(v); out = fill(NaN, n); cs = [0.0; cumsum(v)]
    for i in h+1:n-h
        out[i] = (cs[i+h+1] - cs[i-h]) / (2h + 1)
    end
    return out
end

function study_tolerance(o::StudyOptions; nnoise::Integer = 200)
    F = SGFilter(o.n, o.m; kind = o.kind); m = o.m; g = F.g; S = F.S
    b = nominal_background(); nf = length(b); bmed = median_background()
    inner = m+1:nf-m
    q = noise_pass_fraction(F); S0 = S[m+1]
    d2 = sum(S[k]^2 * (k - m - 1)^2 for k in eachindex(S) if k != m + 1) / q
    Delta = [maximum(abs(b[i+d] - b[i]) for d in -m:m) / b[i] for i in inner]
    b1 = [(b[i+1] - b[i-1]) / 2 for i in inner]
    b2 = [b[i+1] - 2b[i] + b[i-1] for i in inner]
    Deltat = (abs.(b1) .* m .+ 0.5 .* abs.(b2) .* m^2) ./ b[inner]
    sigma0, sigma1 = noise_parameters(Heteroscedastic, o.level, bmed)
    sigma = noise_sigma(sigma0, sigma1, b)
    epsilon = [sum(g[d+m+1]^2 * (sigma[i+d]^2 / sigma[i]^2 - 1) for d in -m:m if d != 0) / q for i in inner]
    epsilont = (sigma1 ./ sigma[inner]) .* (b2 .+ sigma1 .* b1 .^ 2 ./ sigma[inner]) .* d2
    bound = S0 .* (2 .* Delta .+ Delta .^ 2)
    Gb = linear_residual(F, b)
    gb = Dict{Float64,Vector{Float64}}(); dd = Dict{Float64,Vector{Float64}}()
    for lv in o.levels
        s0, s1 = noise_parameters(Heteroscedastic, lv, bmed)
        sigmal = noise_sigma(s0, s1, b)[inner]
        gb[lv] = abs.(Gb) ./ sigmal
        rm = running_mean(Gb ./ sigmal, m)
        sm = zeros(length(rm)); ss = zeros(length(rm))
        for k in 1:nnoise
            x = noise_sigma(s0, s1, b) .* randn(spectrum_rng(o, 960, k), nf)
            u = running_mean(linear_residual(F, x) ./ sigmal, m)
            sm .+= u; ss .+= u .^ 2
        end
        sdv = sqrt.(max.(ss ./ nnoise .- (sm ./ nnoise) .^ 2, 0.0))
        dd[lv] = rm ./ sdv
    end
    res = (; f = collect(frequencies(PAPER_GRID))[inner], Delta = Delta, Delta_taylor = Deltat, eps = epsilon, eps_taylor = epsilont,
           bound, eps_max = running_max(epsilon, m), eps_taylor_max = running_max(epsilont, m), bound_max = running_max(bound, m),
           Gb_sigma = gb, Gb_sigma_max = Dict(k => running_max(v, m) for (k, v) in gb),
           sliding = dd, sliding_max = Dict(k => running_max(replace(v, NaN => 0.0), m) for (k, v) in dd),
           S0, d2, M = kernel_moments(F, F.n % 2 == 0 ? F.n + 2 : F.n + 1), levels = o.levels, level = o.level)
    log_message(@sprintf("tolerance: S0 = %.4f, \u27e8d\u00b2\u27e9 = %.2f, median \u0394 = %.2e, max \u0394 = %.3f, max |\u03b5| = %.2e",
                         S0, d2, median(Delta), maximum(Delta), maximum(abs, epsilon)))
    save_result(result_path(o, "tolerance"), res)
    return res
end
