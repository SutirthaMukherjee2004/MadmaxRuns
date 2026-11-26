Psi(q::Integer, t::Real) = q == 0 ? 1.0 : (t <= 0 ? 0.0 : gamma_inc(q / 2, t / 2)[1])
Psibar(q::Integer, t::Real) = q == 0 ? 0.0 : (t <= 0 ? 1.0 : gamma_inc(q / 2, t / 2)[2])
Phi(x::Real) = 0.5 * erfc(-x / sqrt(2))
Phiinv(p::Real) = -sqrt(2) * erfcinv(2p)

chisq_moment(q::Integer, s::Real) = exp(s * log(2) + loggamma(q / 2 + s) - loggamma(q / 2))

function tilted_tail_moment(q::Integer, t::Real, a::Real, s::Real)
    x = t * (1 + a) / 2
    return (2 / (1 + a))^s * exp(loggamma(q / 2 + s) - loggamma(q / 2)) *
           gamma_inc(q / 2 + s, x)[2] / gamma_inc(q / 2, x)[2]
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
