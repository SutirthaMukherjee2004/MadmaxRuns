Psi(q::Integer, t::Real) = q == 0 ? 1.0 : (t <= 0 ? 0.0 : gamma_inc(q / 2, t / 2)[1])
Psibar(q::Integer, t::Real) = q == 0 ? 0.0 : (t <= 0 ? 1.0 : gamma_inc(q / 2, t / 2)[2])
Phi(x::Real) = 0.5 * erfc(-x / sqrt(2))
Phiinv(p::Real) = -sqrt(2) * erfcinv(2p)

chisq_moment(q::Integer, s::Real) = exp(s * log(2) + loggamma(q / 2 + s) - loggamma(q / 2))

# log of the regularised upper incomplete gamma function Q(a, x). Where Q underflows (deep tail,
# x >> a) the asymptotic series Gamma(a, x) ~ x^(a-1) e^(-x) sum_k (a-1)...(a-k) / x^k is used.
function log_gamma_q(a::Real, x::Real)
    x <= 0 && return 0.0
    Q = gamma_inc(a, x)[2]
    Q > 1e-280 && return log(Q)
    s = 1.0; term = 1.0
    for k in 1:60
        term *= (a - k) / x
        s += term
        abs(term) < 1e-17 * abs(s) && break
    end
    return (a - 1) * log(x) - x - loggamma(a) + log(s)
end
logPsibar(q::Integer, t::Real) = q == 0 ? -Inf : log_gamma_q(q / 2, t / 2)

function tilted_tail_moment(q::Integer, t::Real, a::Real, s::Real)
    x = t * (1 + a) / 2
    return exp(s * log(2 / (1 + a)) + loggamma(q / 2 + s) - loggamma(q / 2) +
               log_gamma_q(q / 2 + s, x) - log_gamma_q(q / 2, x))
end

function chisq_tail_quantile(q::Integer, t::Real, u::Real)
    target = log(u) + log(Psibar(q, t))
    y = max(t, q + 0.0)
    lo, hi = t, max(t, 1.0)
    while log(Psibar(q, hi)) > target
        hi *= 2
    end
    for _ in 1:200
        y = (lo + hi) / 2
        log(Psibar(q, y)) > target ? (lo = y) : (hi = y)
        hi - lo <= 1e-12 * max(1, hi) && break
    end
    return (lo + hi) / 2
end
