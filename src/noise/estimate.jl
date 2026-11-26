struct NoiseFit
    sigma0::Float64
    sigma1::Float64
    bhat::Vector{Float64}
    sigma::Vector{Float64}
    iterations::Int
    converged::Bool
    floor_used::Bool
    sigma0_hist::Vector{Float64}
    sigma1_hist::Vector{Float64}
    estimator::Symbol
    W::Int
end
sigma_profile(f::NoiseFit) = f.sigma

pad_edges(v::AbstractVector, m::Integer) = [fill(v[1], m); v; fill(v[end], m)]

function window_fit(x::AbstractVector, W::Integer, wt::AbstractVector)
    W >= 4 || throw(ArgumentError("W must be at least 4"))
    nw = div(length(x), W)
    nw >= 2 || throw(ArgumentError("need at least two windows"))
    h = (W + 1) / 2
    sw = zeros(nw); bw = zeros(nw)
    @inbounds for j in 1:nw
        o = (j - 1) * W
        S0 = S1 = S2 = S3 = S4 = Y0 = Y1 = Y2 = 0.0
        for k in 1:W
            w = wt[o+k]; t = (k - h) / W; y = x[o+k]; t2 = t * t
            S0 += w; S1 += w * t; S2 += w * t2; S3 += w * t2 * t; S4 += w * t2 * t2
            Y0 += w * y; Y1 += w * t * y; Y2 += w * t2 * y
        end
        A = [S0 S1 S2; S1 S2 S3; S2 S3 S4]
        c0, c1, c2 = A \ [Y0, Y1, Y2]
        r2 = 0.0; fb = 0.0
        for k in 1:W
            t = (k - h) / W; fk = c0 + c1 * t + c2 * t * t
            r2 += (x[o+k] - fk)^2; fb += fk
        end
        sw[j] = sqrt(r2 / (W - 3)); bw[j] = fb / W
    end
    bm = mean(bw); sm = mean(sw)
    s1 = sum((bw .- bm) .* (sw .- sm)) / sum(abs2, bw .- bm)
    return sm - s1 * bm, s1, sw, bw
end
function lsq_noise_fit(r::AbstractVector, bhat::AbstractVector, q::Real;
                       init = (sqrt(mean(abs2, r) / q), 0.0), maxit::Integer = 200, tol::Real = 1e-13)
    p = collect(Float64, init); y = r .^ 2
    cost(p) = sum(abs2, y .- q .* (p[1] .+ p[2] .* bhat) .^ 2)
    lambda = 1e-3; c = cost(p)
    bsc = max(mean(abs, bhat), eps())
    for _ in 1:maxit
        s = p[1] .+ p[2] .* bhat
        res = y .- q .* s .^ 2
        J1 = 2q .* s; J2 = 2q .* s .* bhat
        A = [sum(J1 .* J1) sum(J1 .* J2); sum(J1 .* J2) sum(J2 .* J2)]
        g = [sum(J1 .* res), sum(J2 .* res)]
        improved = false
        for _ in 1:30
            Delta = (A + lambda * Diagonal(diag(A))) \ g
            pn = p .+ Delta; cn = cost(pn)
            if cn < c
                done = abs(Delta[1]) + abs(Delta[2]) * bsc <= tol * (abs(pn[1]) + abs(pn[2]) * bsc)
                p, c = pn, cn; lambda = max(lambda / 10, 1e-12); improved = true
                done && return p[1], p[2]
                break
            end
            lambda *= 10
        end
        improved || break
    end
    return p[1], p[2]
end

function estimate_noise(F::SGFilter, cache::WSGCache, x::AbstractVector; W::Integer = 96,
                        estimator::Symbol = :window, tol::Real = 1e-10, maxit::Integer = 60,
                        floor_frac::Real = 0.1)
    estimator in (:window, :lsq) || throw(ArgumentError("estimator must be :window or :lsq"))
    m = F.m; nf = length(x); N = nf - 2m
    q = noise_pass_fraction(F)
    xin = view(x, m+1:nf-m)
    bhat = smooth_unweighted(F, x)
    wt = ones(nf)
    fitpars(bh, wt) = estimator == :window ? window_fit(x, W, wt)[1:2] :
                      lsq_noise_fit(xin .- bh, bh, q)
    s0, s1 = fitpars(bhat, wt)
    h0 = [s0]; h1 = [s1]
    bsc = median(abs.(bhat))
    conv = false; flused = false; nit = 0
    sig = s0 .+ s1 .* bhat
    for k in 1:maxit
        sf = pad_edges(sig, m)
        fl = floor_frac * mean(abs, sf)
        @inbounds for i in 1:nf
            v = sf[i]
            v < fl && (v = fl; flused = true)
            wt[i] = 1 / v^2
        end
        bhat = weighted_smooth(cache, x, wt)
        n0, n1 = estimator == :window ? window_fit(x, W, wt)[1:2] : lsq_noise_fit(xin .- bhat, bhat, q; init = (s0, s1))
        sc = abs(s0) + abs(s1) * bsc
        conv = abs(n0 - s0) <= tol * sc && abs(n1 - s1) * bsc <= tol * sc
        s0, s1 = n0, n1; push!(h0, s0); push!(h1, s1); nit = k
        sig = s0 .+ s1 .* bhat
        conv && break
    end
    return NoiseFit(s0, s1, bhat, sig, nit, conv, flused, h0, h1, estimator, W)
end
