squares_T(z::AbstractVector) = Float64(RunStatistics.t_obs(collect(Float64, z), 0.0, 1.0)[1])

function squares_T_run(z::AbstractVector)
    T = 0.0; s = 0.0; i0 = 0; best = (0, 0)
    @inbounds for i in eachindex(z)
        if z[i] >= 0
            s == 0 && (i0 = i)
            s += z[i]^2
            s > T && (T = s; best = (i0, i))
        else
            s = 0.0
        end
    end
    return T, best[1], best[2]
end

function pbc_block_logcdfs(t::Real, L::Integer)
    weights = [exp2(-k) for k in 0:L]
    powers = [gamma_inc(k / 2, t / 2) for k in 1:L]
    cdf_positive = zeros(L + 1)
    survival = zeros(L + 1)
    logs = zeros(L)
    for n in 1:L
        lo, hi = powers[n]
        lower = 0.5 * cdf_positive[n] + weights[n+1] * lo
        upper = 0.5 * survival[n] + weights[n+1] * hi
        @inbounds for k in 1:n-1
            lo, hi = powers[k]
            w = weights[k+2]
            lower += w * lo * (cdf_positive[n-k] + weights[n-k])
            upper += w * (hi + lo * survival[n-k])
        end
        cdf_positive[n+1] = lower
        survival[n+1] = upper
        normalizer = 1 - weights[n+1]
        logs[n] = lower <= 0.5 * normalizer ? log(lower / normalizer) : log1p(-upper / normalizer)
    end
    return logs
end

function pbc_logcdf(t::Real, L::Integer; nblock::Integer = 100, epsp::Real = 1e-7)
    L > 0 || throw(ArgumentError("spectrum length must be positive"))
    nblock > 0 || throw(ArgumentError("block length must be positive"))
    isfinite(t) || throw(DomainError(t, "SQUARES statistic must be finite"))
    t <= 0 && return -Inf
    L <= 2nblock && return pbc_block_logcdfs(t, L)[end]
    return pbc_logcdf_blocks(t, L, nblock)
end
pbc_cdf(t::Real, L::Integer; kw...) = exp(pbc_logcdf(t, L; kw...))

function checked_pvalue(p::Real, t::Real)
    0 < p <= 1 || throw(DomainError(p, "invalid or unrepresentable p-value at T=$t"))
    return p
end

pbc_pvalue(t::Real, L::Integer; kw...) = checked_pvalue(-expm1(pbc_logcdf(t, L; kw...)), t)

function pbc_cdf_recursion(t::Real, L::Integer; qmax::Integer = 400)
    qm = min(L, qmax)
    w = [0.5^(k + 1) * Psi(k, t) for k in 0:qm]
    a = zeros(L + 1); a[1] = 1.0
    @inbounds for n in 1:L
        s = 0.0
        for k in 0:min(n - 1, qm)
            s += a[n-k] * w[k+1]
        end
        a[n+1] = s
    end
    P = 0.0
    @inbounds for k in 0:qm
        P += a[L-k+1] * 0.5^k * Psi(k, t)
    end
    return P
end

function pbc_logcdf_blocks(t::Real, L::Integer, nblock::Integer)
    blocks = pbc_block_logcdfs(t, 2nblock)
    F1, F2 = blocks[nblock], blocks[2nblock]
    n = L / nblock
    return F1 == -Inf ? -Inf : F1 + (n - 1) * (F2 - F1)
end

pbc_cdf_blocks(t::Real, L::Integer; nblock::Integer = 100) = t <= 0 ? 0.0 : exp(pbc_logcdf_blocks(t, L, nblock))
