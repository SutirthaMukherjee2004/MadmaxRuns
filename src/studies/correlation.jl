const CORR_TAUS = [1e-5, 1e-4, 1e-3, 0.01, 0.05, 0.2]
const CORR_BLOCK = 160
const TAUSCAN_LAGS = 250

corr_block(Ng::Integer, L::Integer = CORR_BLOCK) = (div(Ng, 2) - div(L, 2) + 1):(div(Ng, 2) - div(L, 2) + L)

function autocorrelation(z::AbstractVector, L::Integer)
    n = length(z); P = nextprod([2, 3, 5], 2n)
    Z = rfft([z; zeros(P - n)])
    a = irfft(abs2.(Z), P)
    return a[2:L+1] ./ a[1]
end

accumulate_term(a, b) = a + b
accumulate_term(a::AbstractArray, b::AbstractArray) = (a .+= b)

function tsum(f, n::Integer)
    nt = min(Threads.nthreads(), n)
    chunks = [i:nt:n for i in 1:nt]
    parts = Vector{Any}(undef, nt)
    Threads.@threads for c in 1:nt
        acc = nothing
        for i in chunks[c]
            v = f(i)
            acc = acc === nothing ? deepcopy(v) : map(accumulate_term, acc, v)
        end
        parts[c] = acc
    end
    return reduce((a, b) -> map(accumulate_term, a, b), filter(!isnothing, parts))
end

cov2corr(C::AbstractMatrix) = (d = sqrt.(diag(C)); C ./ (d * d'))

function mean_lag(C::AbstractMatrix, L::Integer = 10)
    s = 0.0; k = 0
    for i in axes(C, 1), j in axes(C, 2)
        1 <= abs(i - j) <= L && (s += C[i, j]; k += 1)
    end
    return s / k
end

function iid_reference(S::FilterSetup, n::Integer; seed::Integer, lags::Integer = TAUSCAN_LAGS,
                       block = corr_block(S.Ng))
    acc = tsum(n) do i
        Y = rfft(randn(stream_rng(seed, i), S.N))
        Y[.!S.keep] .= 0
        v = irfft(Y, S.N)[S.rows] ./ sqrt(S.kappa)
        b = v[block]
        (b * b', autocorrelation(v, lags))
    end
    return cov2corr(acc[1]), acc[2] ./ n
end

function corr_terms(z, block)
    b = z[block]
    return (b * b', autocorrelation(z, TAUSCAN_LAGS))
end

function study_correlation(o::StudyOptions)
    tmin, t0h, t0x = thresholds(o)
    taus = sort(unique([CORR_TAUS; tmin; t0h; t0x]))
    cfg = base_config(o)
    setups = tau_setups(cfg, taus; with_pk = false, with_boundary = false)
    S1 = setups[1]; blk = corr_block(S1.Ng)
    refs = Dict{Float64,Any}()
    for (k, S) in enumerate(setups)
        refs[taus[k]] = iid_reference(S, o.nnull; seed = o.seed + 31 + k)
    end
    cases = [(model, lv, true) for model in o.models for lv in o.levels]
    push!(cases, (Homoscedastic, o.level, false))
    out = Dict{Tuple{String,Float64,Bool},Any}()
    for (model, lv, bon) in cases
        gen = data_generator(o, cfg, case_id(model, lv; background_on = bon); model, level = lv, background_on = bon)
        set = fit_set(S1, "corr", gen, o.nspec)
        terms = tsum(length(set)) do i
            sigma = set.sigma[i]
            r = linear_residual(S1.F, set.gen(i))
            R = rfft(r)
            raw = r[S1.rows] ./ (sigma[S1.rows] .* sqrt(S1.q))
            whitened = (corr_terms(tau_whitened(S, sigma, R), blk) for S in setups)
            (corr_terms(raw, blk)..., Iterators.flatten(whitened)...)
        end
        Craw = cov2corr(terms[1])
        pertau = Dict{Float64,Any}()
        for (k, S) in enumerate(setups)
            C = cov2corr(terms[2k+1]); rho = terms[2k+2] ./ length(set)
            Cref, rhoref = refs[taus[k]]
            D = C .- Cref; D[diagind(D)] .= 0
            pertau[taus[k]] = (; corr = C, diff = D, mean_lag = mean_lag(C), mean_diff = mean_lag(D),
                              delta_rms = sqrt(mean(abs2, rho .- rhoref)), rho = rho)
        end
        out[(model_tag(model), lv, bon)] = (; raw = Craw, raw_mean_lag = mean_lag(Craw), rho_raw = terms[2] ./ length(set), taus = pertau)
        log_message(@sprintf("correlation %s %.2f%% b=%s: raw mean lag 1-10 = %.4f; \u0394_rms at \u03c4_min = %.4f",
                             model_tag(model), 100lv, bon ? "on" : "0", mean_lag(Craw), pertau[tmin].delta_rms))
    end
    res = (; taus, figure_taus = CORR_TAUS, tau_min = tmin, tau0 = (homo = t0h, hetero = t0x), block = blk,
           cases = out, refs, levels = o.levels, level = o.level, models = model_tag.(o.models), nspec = o.nspec, nnull = o.nnull)
    save_result(result_path(o, "correlation"), res)
    return res
end

function study_tauscan(o::StudyOptions)
    taus = tau_grid(o); cfg = base_config(o)
    setups = tau_setups(cfg, taus; with_pk = false, with_boundary = false)
    S1 = setups[1]
    rmsref = zeros(length(taus))
    for (k, S) in enumerate(setups)
        _, rho = iid_reference(S, o.nnull; seed = o.seed + 500 + k)
        rmsref[k] = sqrt(mean(abs2, rho))
    end
    out = Dict{Tuple{String,Float64},Vector{Float64}}()
    for model in o.models, lv in o.levels
        set = fit_set(S1, "tauscan", data_generator(o, cfg, case_id(model, lv); model, level = lv), o.nspec)
        rhos = tsum(length(set)) do i
            R = rfft(linear_residual(S1.F, set.gen(i)))
            sigma = set.sigma[i]
            Tuple(autocorrelation(tau_whitened(S, sigma, R), TAUSCAN_LAGS) for S in setups)
        end
        rms = [sqrt(mean(abs2, rho ./ length(set))) for rho in rhos]
        out[(model_tag(model), lv)] = rms
        log_message(@sprintf("tauscan %s %.2f%%: minimum rms %.5f at \u03c4 = %.2e", model_tag(model), 100lv,
                             minimum(rms), taus[argmin(rms)]))
    end
    mrms = mean(values(out))
    tmin = parabolic_minimum(log10.(taus), mrms)
    log_message(@sprintf("tauscan: \u03c4_min = %.2e", tmin))
    res = (; taus, rms = out, rms_ref = rmsref, tau_min = tmin, levels = o.levels, models = model_tag.(o.models),
           nspec = o.nspec, nnull = o.nnull)
    save_result(result_path(o, "tauscan"), res)
    return res
end

function parabolic_minimum(x::AbstractVector, y::AbstractVector)
    k = argmin(y)
    (k == 1 || k == length(y)) && return exp10(x[k])
    x1, x2, x3 = x[k-1], x[k], x[k+1]; y1, y2, y3 = y[k-1], y[k], y[k+1]
    d = (x1 - x2) * (x1 - x3) * (x2 - x3)
    a = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / d
    b = (x3^2 * (y1 - y2) + x2^2 * (y3 - y1) + x1^2 * (y2 - y3)) / d
    a <= 0 && return exp10(x2)
    return exp10(clamp(-b / (2a), x1, x3))
end
