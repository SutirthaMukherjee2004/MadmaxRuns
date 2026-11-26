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

function pbc_cdf(t::Real, L::Integer; nblock::Integer = 100, epsp::Real = 1e-7)
    t <= 0 && return 0.0
    L <= 2nblock && return Float64(RunStatistics.squares_cdf(Float64(t), L))
    for e in (epsp, 1e-10, 1e-15)
        try
            return 1 - Float64(RunStatistics.squares_pvalue_approx(Float64(t), Float64[nblock, L / nblock], e))
        catch err
            (err isa AssertionError || err isa DomainError) || rethrow()
        end
    end
    return 1 - 1e-15
end
pbc_pvalue(t::Real, L::Integer; kw...) = 1 - pbc_cdf(t, L; kw...)

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

function pbc_cdf_blocks(t::Real, L::Integer; nblock::Integer = 100)
    F1 = pbc_cdf_recursion(t, nblock); F2 = pbc_cdf_recursion(t, 2nblock)
    n = L / nblock
    return F1^n / (F1^2 / F2)^(n - 1)
end
