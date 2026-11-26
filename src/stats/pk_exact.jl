function Eq_exact(gamma::AbstractVector, q::Integer, t::Real; M::Integer = 200_000,
                  rng::AbstractRNG = Random.default_rng())
    R = correlation_matrix(gamma, q); W = inv(R) - I
    d = det(R)^-0.5
    w = zeros(q + 2); g = zeros(q)
    s = 0.0; s2 = 0.0
    for _ in 1:M
        Y = chisq_tail_quantile(q, t, rand(rng))
        randn!(rng, g); g .= abs.(g); g ./= norm(g)
        w[1] = -abs(randn(rng)); w[end] = -abs(randn(rng))
        @views w[2:q+1] .= sqrt(Y) .* g
        x = exp(-0.5 * dot(w, W * w))
        s += x; s2 += x^2
    end
    m = s / M
    return d * m, d * sqrt(max(s2 / M - m^2, 0.0) / M)
end
