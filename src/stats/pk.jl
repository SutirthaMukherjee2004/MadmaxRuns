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
    run = (1 + a)^(-q / 2) * Psibar(q, t * (1 + a)) / Psibar(q, t)
    nb(eta) = 2 * exp(eta^2 / 2) * Phi(eta)
    E = d.detfac * run * nb(etam) * nb(etap)
    return isfinite(E) ? E : 0.0
end

function Eq_second(d::EqData, t::Real)
    q, W, mo = d.q, d.W, d.mom
    a2 = mo.Eu - t * mo.Vuu / 2
    1 + a2 <= 0 && return 0.0
    run = (1 + a2)^(-q / 2) * Psibar(q, t * (1 + a2)) / Psibar(q, t) * exp(-t^2 * mo.Vuu / 8)
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
    return isfinite(E) ? E : 0.0
end

struct PKModel
    N::Int
    Ng::Int
    K::Vector{Int}
    kappa::Float64
    gamma::Vector{Float64}
    eq::Vector{EqData}
    order::Int
end
function PKModel(N::Integer, Ng::Integer, K::AbstractVector{<:Integer}; qmax::Integer = 60, order::Integer = 2)
    order in (1, 2) || throw(ArgumentError("order must be 1 or 2"))
    gamma, kappa = gamma_lags(N, K, qmax + 1)
    eq = Vector{EqData}(undef, qmax)
    Threads.@threads for q in 1:qmax
        eq[q] = EqData(gamma, q)
    end
    return PKModel(N, Ng, collect(K), kappa, gamma, eq, order)
end

function pk_cdf(M::PKModel, t::Real; Eqf = nothing, P0 = nothing)
    t <= 0 && return 0.0
    S = 0.0
    for d in M.eq
        pi0 = 0.5^(d.q + 2) * Psibar(d.q, t)
        pi0 < 1e-300 && continue
        E = Eqf !== nothing ? Eqf(d.q, t) : (M.order == 2 ? Eq_second(d, t) : Eq_first(d, t))
        S += pi0 * (E - 1)
    end
    pb = P0 === nothing ? pbc_cdf(t, M.Ng) : P0
    return pb * exp(-M.Ng * S)
end
pk_pvalue(M::PKModel, t::Real; kw...) = 1 - pk_cdf(M, t; kw...)

struct PKTable
    t::Vector{Float64}
    logS::Vector{Float64}
end
function PKTable(cdf::Function; tmax::Real = 150.0, dt::Real = 0.05)
    t = collect(0.0:dt:tmax)
    S = zeros(length(t))
    Threads.@threads for i in eachindex(t)
        S[i] = 1 - cdf(t[i])
    end
    return PKTable(t, log.(clamp.(S, 1e-300, 1.0)))
end
pk_table(M::PKModel; kw...) = PKTable(t -> pk_cdf(M, t); kw...)
pbc_table(L::Integer; kw...) = PKTable(t -> pbc_cdf(t, L); kw...)

function pvalue(tab::PKTable, T::Real)
    T <= tab.t[1] && return exp(tab.logS[1])
    if T >= tab.t[end]
        k = length(tab.t); sl = (tab.logS[k] - tab.logS[k-1]) / (tab.t[k] - tab.t[k-1])
        return exp(tab.logS[k] + sl * (T - tab.t[k]))
    end
    k = searchsortedlast(tab.t, T); a = (T - tab.t[k]) / (tab.t[k+1] - tab.t[k])
    return exp((1 - a) * tab.logS[k] + a * tab.logS[k+1])
end

function pk_threshold(cdf::Function, alpha::Real; lo::Real = 1.0, hi::Real = 150.0, tol::Real = 1e-7)
    while hi - lo > tol
        mid = (lo + hi) / 2
        cdf(mid) < 1 - alpha ? (lo = mid) : (hi = mid)
    end
    return (lo + hi) / 2
end
