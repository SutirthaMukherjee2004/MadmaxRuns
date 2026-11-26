function gauss_legendre(n::Integer)
    beta = [k / sqrt(4k^2 - 1) for k in 1:n-1]
    E = eigen(SymTridiagonal(zeros(n), beta))
    return E.values, 2 .* E.vectors[1, :] .^ 2
end
const GL_NODES, GL_WEIGHTS = gauss_legendre(48)

function Phi2(h::Real, k::Real, rho::Real)
    abs(rho) <= 1 || throw(DomainError(rho, "correlation must lie in [-1, 1]"))
    base = Phi(h) * Phi(k)
    rho == 0 && return base
    if isinf(h) || isinf(k)
        (h == -Inf || k == -Inf) && return 0.0
        return Phi(min(h, k))
    end
    thetamax = asin(rho); npan = 4; s = 0.0
    for p in 0:npan-1
        a = thetamax * p / npan; b = thetamax * (p + 1) / npan
        hw = (b - a) / 2; mid = (a + b) / 2
        for (x, w) in zip(GL_NODES, GL_WEIGHTS)
            theta = mid + hw * x
            c2 = cos(theta)^2
            s += w * hw * exp(-(h^2 - 2h * k * sin(theta) + k^2) / (2c2))
        end
    end
    return base + s / (2 * pi)
end
