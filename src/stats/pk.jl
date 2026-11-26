function gamma_lags(N::Integer, K::AbstractVector{<:Integer}, ellmax::Integer)
    kappa = 1 - sum((mode_weight(f, N) for f in K); init = 0) / N
    gamma = [-sum((mode_weight(f, N) * cos((2 * pi) * f * ell / N) for f in K); init = 0.0) / (N * kappa) for ell in 1:ellmax]
    return gamma, kappa
end

function correlation_matrix(gamma::AbstractVector, q::Integer)
    n = q + 2
    R = Matrix{Float64}(I, n, n)
    for a in 1:n, b in 1:n
        a != b && (R[a, b] = gamma[abs(a - b)])
    end
    return R
end

const HN_MEAN = sqrt(2 / pi)
const HN_VAR = 1 - 2 / pi
const HN_MU3 = HN_MEAN * (4 / pi - 1)
const HN_MU4 = 3 - 2 * (2 / pi) - 3 * (2 / pi)^2

function direction_moments(A::AbstractMatrix, b::AbstractVector, c::AbstractVector)
    q = size(A, 1)
    mu1, sigma2, mu3, mu4 = HN_MEAN, HN_VAR, HN_MU3, HN_MU4
    one = ones(q); dA = diag(A); Amu = A * (mu1 .* one)
    EQ = sigma2 * tr(A) + mu1^2 * sum(A)
    VQ = (mu4 - 3 * sigma2^2) * sum(abs2, dA) + (2 * sigma2^2) * sum(abs2, A) + (4 * sigma2) * dot(Amu, Amu) + (4 * mu3) * dot(Amu, dA)
    EQ2 = VQ + EQ^2
    Lb = mu1 * sum(b); Lc = mu1 * sum(c)
    EQLb = EQ * Lb + (2 * sigma2) * dot(Amu, b) + mu3 * dot(dA, b)
    EQLc = EQ * Lc + (2 * sigma2) * dot(Amu, c) + mu3 * dot(dA, c)
    ELbb = sigma2 * dot(b, b) + Lb^2; ELbc = sigma2 * dot(b, c) + Lb * Lc; ELcc = sigma2 * dot(c, c) + Lc^2
    Y_half, Y1, Y_three_halves, Y2 = chisq_moment(q, 0.5), chisq_moment(q, 1), chisq_moment(q, 1.5), chisq_moment(q, 2)
    Eu = EQ / Y1; Eb = Lb / Y_half; Ec = Lc / Y_half
    return (Eu = Eu, Eb = Eb, Ec = Ec, Vuu = EQ2 / Y2 - Eu^2,
            Cub = EQLb / Y_three_halves - Eu * Eb, Cuc = EQLc / Y_three_halves - Eu * Ec,
            Cbb = ELbb / Y1 - Eb^2, Cbc = ELbc / Y1 - Eb * Ec, Ccc = ELcc / Y1 - Ec^2)
end

struct EqData
    q::Int
    detfac::Float64
    W::Matrix{Float64}
    a1::Float64
    c1::Float64
    mom::NamedTuple
end
function EqData(gamma::AbstractVector, q::Integer)
    R = correlation_matrix(gamma, q)
    W = inv(R) - I
    Wvv = W[2:q+1, 2:q+1]
    a1 = tr(Wvv) / q + 2 / (pi * q) * (sum(Wvv) - tr(Wvv))
    c1 = exp(loggamma(q / 2) - loggamma((q + 1) / 2)) / sqrt(pi)
    mom = direction_moments(Wvv, W[1, 2:q+1], W[q+2, 2:q+1])
    return EqData(q, det(R)^-0.5, W, a1, c1, mom)
end

function Eq_first(d::EqData, t::Real)
    q, W = d.q, d.W
    a = d.a1
    1 + a <= 0 && return 0.0
    etam = sqrt(t) * d.c1 * sum(view(W, 1, 2:q+1))
    etap = sqrt(t) * d.c1 * sum(view(W, q + 2, 2:q+1))
    run = exp(-q / 2 * log1p(a) + logPsibar(q, t * (1 + a)) - logPsibar(q, t))
    nb(eta) = 2 * exp(eta^2 / 2) * Phi(eta)
    E = d.detfac * run * nb(etam) * nb(etap)
    return E
end

function Eq_second(d::EqData, t::Real)
    q, W, mo = d.q, d.W, d.mom
    a2 = mo.Eu - t * mo.Vuu / 2
    1 + a2 <= 0 && return 0.0
    run = exp(-q / 2 * log1p(a2) + logPsibar(q, t * (1 + a2)) - logPsibar(q, t) - t^2 * mo.Vuu / 8)
    y_half = tilted_tail_moment(q, t, a2, 0.5)
    y1 = tilted_tail_moment(q, t, a2, 1.0)
    y_three_halves = tilted_tail_moment(q, t, a2, 1.5)
    eta = [y_half * mo.Eb - 0.5 * y_three_halves * mo.Cub, y_half * mo.Ec - 0.5 * y_three_halves * mo.Cuc]
    C = y1 .* [mo.Cbb mo.Cbc; mo.Cbc mo.Ccc] .- [W[1, 1] W[1, q+2]; W[q+2, 1] W[q+2, q+2]]
    P = Matrix{Float64}(I, 2, 2) .- C
    isposdef(Symmetric(P)) || return 0.0
    Pi = inv(P); mu = Pi * eta
    s = sqrt.(diag(Pi)); rho = clamp(Pi[1, 2] / (s[1] * s[2]), -1.0, 1.0)
    nb = 4 / sqrt(det(P)) * exp(0.5 * dot(eta, mu)) * Phi2(mu[1] / s[1], mu[2] / s[2], rho)
    E = d.detfac * run * nb
    return E
end

# Appendix G: P_K(T < t) = P_BC(T < t | N_g) exp[-N_g sum_q pi0_q (E_q - 1)], pi0_q = 2^-(q+2) Psibar_q(t).
#
# Two numerical rules complete Eq. (G.4); both are applied identically in every caller.
#
# 1. Long runs. Where the closed form for E_q is undefined (P not positive definite, or 1 + a'_q -> 0)
#    and returns 0, NaN or Inf, E_q = 0, as stated below Eq. (G.5). This is accepted only when the
#    dropped run lengths carry negligible weight, sum_dropped pi0_q <= PK_DROP_REL * sum_q pi0_q E_q.
#    Since E_q <= 1 for long runs, they then change p_K by at most 0.01 %. Otherwise it is an error.
#
# 2. Deep tail. Eq. (G.4) is validated against Monte Carlo down to p = 1e-4 (App. G). Below that its
#    CDF form subtracts nearly equal numbers and the second-order E_q leaves its range of validity
#    (1 + a'_q -> 0 as t grows). Beyond the switch point t*, where Eq. (G.4) gives p = PK_TAIL_SWITCH
#    (or the last T at which rule 1 holds, if that comes first), the correction factor is frozen:
#        p_K(t) = p_BC(t) * p_K(t*) / p_BC(t*),   t > t*.
#    P_BC is exact (Eq. G.1). The factor p_K/p_BC decreases with t (removing the lowest modes makes
#    long runs rarer, Eq < 1), so the frozen factor never claims more significance than Eq. (G.4) would.
const PK_TAIL_SWITCH = 1e-4
const PK_DROP_REL = 1e-4
const PK_TAIL_STEP = 0.25

struct PKModel
    N::Int
    Ng::Int
    K::Vector{Int}
    kappa::Float64
    gamma::Vector{Float64}
    eq::Vector{EqData}
    order::Int
    tail_t::Float64
    tail_logratio::Float64
end
function PKModel(N::Integer, Ng::Integer, K::AbstractVector{<:Integer}; qmax::Integer = 60, order::Integer = 2)
    order in (1, 2) || throw(ArgumentError("order must be 1 or 2"))
    gamma, kappa = gamma_lags(N, K, qmax + 1)
    eq = Vector{EqData}(undef, qmax)
    Threads.@threads for q in 1:qmax
        eq[q] = EqData(gamma, q)
    end
    M = PKModel(N, Ng, collect(K), kappa, gamma, eq, order, Inf, 0.0)
    tail_t, tail_logratio = pk_tail_match(M)
    return PKModel(N, Ng, collect(K), kappa, gamma, eq, order, tail_t, tail_logratio)
end

function logsumexp_acc(a::Real, b::Real)
    a == -Inf && return b
    b == -Inf && return a
    m = max(a, b)
    return m + log1p(exp(-abs(a - b)))
end

# S = sum_q pi0_q (E_q - 1) with rule 1. Returns nothing (non-strict) or throws (strict) when a
# non-negligible E_q is invalid.
function pk_terms(M::PKModel, t::Real, Eqf, strict::Bool)
    S = 0.0; lsE = -Inf
    dropped_lpi = -Inf; dropped_q = 0
    for d in M.eq
        lpi = -(d.q + 2) * log(2.0) + logPsibar(d.q, t)
        lpi == -Inf && continue
        E = Eqf !== nothing ? Eqf(d.q, t) : (M.order == 2 ? Eq_second(d, t) : Eq_first(d, t))
        if !isfinite(E) || E < 0
            dropped_lpi = logsumexp_acc(dropped_lpi, lpi); dropped_q = max(dropped_q, d.q)
            S -= exp(lpi)
            continue
        end
        S += exp(lpi) * (E - 1)
        E > 0 && (lsE = logsumexp_acc(lsE, lpi + log(E)))
    end
    if dropped_lpi > -Inf && !(dropped_lpi <= log(PK_DROP_REL) + lsE)
        strict && throw(DomainError(t, "invalid P_K run correction E_q (up to q = $dropped_q) with non-negligible weight at T = $t"))
        return nothing
    end
    return S
end

function pk_logcdf(M::PKModel, t::Real, Eqf, P0, strict::Bool)
    isfinite(t) || throw(DomainError(t, "SQUARES statistic must be finite"))
    t <= 0 && return -Inf
    S = pk_terms(M, t, Eqf, strict)
    S === nothing && return NaN
    P0 === nothing || 0 <= P0 <= 1 || throw(DomainError(P0, "P_BC must be a probability"))
    baseline = P0 === nothing ? pbc_logcdf(t, M.Ng) : log(P0)
    correction = M.Ng * S
    result = baseline - correction
    if !isfinite(correction) || !(result <= 0)
        strict && throw(DomainError(result, "invalid log P_K CDF at T=$t"))
        return NaN
    end
    return result
end

# log P_BC(T >= t | L). Where p_BC underflows the exact expression (T >~ 1400), the leading Poisson
# term L sum_q 2^-(q+2) Psibar_q(t) is used; it agrees with the exact tail to 1e-3 from T ~ 40 on.
function pbc_logsurvival(t::Real, L::Integer)
    l = pbc_logcdf(t, L)
    l < 0 && return log(-expm1(l))
    ls = -Inf
    for q in 1:min(L, 400)
        ls = logsumexp_acc(ls, -(q + 2) * log(2.0) + logPsibar(q, t))
    end
    return log(L) + ls
end

function pk_tail_match(M::PKModel)
    target = log(PK_TAIL_SWITCH)
    logp(t) = (l = pk_logcdf(M, t, nothing, nothing, false); isfinite(l) && l < 0 ? log(-expm1(l)) : NaN)
    t = PK_TAIL_STEP; prev = t; lp = logp(t)
    isfinite(lp) || error("P_K is not defined at T = $t")
    while lp > target
        nt = t + PK_TAIL_STEP
        nl = logp(nt)
        isfinite(nl) || break                       # rule 1 fails first: freeze at the last valid T
        prev, t, lp = t, nt, nl
        t > 1e4 && error("P_K survival does not reach $(PK_TAIL_SWITCH)")
    end
    if lp <= target                                 # refine the crossing inside [prev, t]
        lo, hi = prev, t
        for _ in 1:60
            mid = (lo + hi) / 2
            lm = logp(mid)
            (isfinite(lm) && lm > target) ? (lo = mid) : (hi = mid)
            hi - lo <= 1e-10 * hi && break
        end
        t = lo; lp = logp(t)
    end
    return t, lp - pbc_logsurvival(t, M.Ng)
end

# log P_K(T >= t): Eq. (G.4) up to the switch point, the frozen correction factor beyond it.
function pk_logsurvival(M::PKModel, t::Real)
    isfinite(t) || throw(DomainError(t, "SQUARES statistic must be finite"))
    t <= 0 && return 0.0
    if t <= M.tail_t
        l = pk_logcdf(M, t, nothing, nothing, true)
        l < 0 || throw(DomainError(l, "invalid log P_K CDF at T=$t"))
        return log(-expm1(l))
    end
    return min(0.0, pbc_logsurvival(t, M.Ng) + M.tail_logratio)
end

pk_cdf(M::PKModel, t::Real; Eqf = nothing, P0 = nothing) =
    Eqf === nothing && P0 === nothing ? -expm1(pk_logsurvival(M, t)) : exp(pk_logcdf(M, t, Eqf, P0, true))
# p values below the smallest normal Float64 are reported as floatmin(Float64); pk_log10pvalue keeps them.
pk_pvalue(M::PKModel, t::Real; Eqf = nothing, P0 = nothing) =
    Eqf === nothing && P0 === nothing ? checked_pvalue(max(exp(pk_logsurvival(M, t)), floatmin(Float64)), t) :
    checked_pvalue(-expm1(pk_logcdf(M, t, Eqf, P0, true)), t)
pk_log10pvalue(M::PKModel, t::Real) = pk_logsurvival(M, t) / log(10)

struct PKTable
    t::Vector{Float64}
    logS::Vector{Float64}
end
function survival_table(survival::Function; tmax::Real = 150.0, dt::Real = 0.05, allow_invalid::Bool = false)
    isfinite(dt) && isfinite(tmax) && 0 < dt <= tmax ||
        throw(ArgumentError("p-value table requires finite 0 < dt <= tmax"))
    t = collect(0.0:dt:tmax)
    logS = zeros(length(t))
    Threads.@threads for i in eachindex(t)
        p = survival(t[i])
        logS[i] = allow_invalid && !(0 < p <= 1) ? NaN : log(checked_pvalue(p, t[i]))
    end
    return PKTable(t, logS)
end
PKTable(cdf::Function; kw...) = survival_table(t -> 1 - cdf(t); kw...)
pk_table(M::PKModel; kw...) = survival_table(t -> exp(pk_logsurvival(M, t)); kw...)
pbc_table(L::Integer; kw...) = survival_table(t -> pbc_pvalue(t, L); kw...)

function pvalue(tab::PKTable, T::Real)
    isfinite(T) || throw(DomainError(T, "SQUARES statistic must be finite"))
    T <= tab.t[1] && return checked_pvalue(exp(tab.logS[1]), T)
    T <= last(tab.t) || throw(DomainError(T, "SQUARES statistic exceeds the lookup table range"))
    k = searchsortedlast(tab.t, T)
    T == tab.t[k] && return checked_pvalue(exp(tab.logS[k]), T)
    a = (T - tab.t[k]) / (tab.t[k+1] - tab.t[k])
    logS = (1 - a) * tab.logS[k] + a * tab.logS[k+1]
    return checked_pvalue(exp(logS), T)
end

function tabulated_pvalue(evaluate, tab::Union{Nothing,PKTable}, T::Real)
    return tab === nothing || T > last(tab.t) ? evaluate(T) : pvalue(tab, T)
end

function decision_threshold(tab::PKTable, alpha::Real)
    0 < alpha < 1 || throw(ArgumentError("significance must be in (0, 1)"))
    target = log(alpha)
    for k in eachindex(tab.t)
        logp = tab.logS[k]
        isfinite(logp) || error("invalid p-value table before reaching significance $alpha")
        logp > target && continue
        k == 1 && return tab.t[k]
        before = tab.logS[k-1]
        weight = (target - before) / (logp - before)
        return tab.t[k-1] + weight * (tab.t[k] - tab.t[k-1])
    end
    error("p-value table does not reach significance $alpha")
end

function pk_threshold(cdf::Function, alpha::Real; lo::Real = 1.0, hi::Real = 150.0, tol::Real = 1e-7)
    while hi - lo > tol
        mid = (lo + hi) / 2
        cdf(mid) < 1 - alpha ? (lo = mid) : (hi = mid)
    end
    return (lo + hi) / 2
end
