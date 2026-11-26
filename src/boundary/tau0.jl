function boundary_significance(F::SGFilter, Dspk::AbstractMatrix, G::AbstractVector, sigmafull::AbstractVector,
                               tau::Real; alpha::Real = 0.05, drop_nyquist::Bool = true)
    m = F.m; nfull = length(sigmafull); N = nfull - 2m
    keep = retained_modes(G, tau; drop_nyquist, N); kappa = retained_fraction(keep, N)
    rows = tested_bins(N, m); Ng = length(rows)
    J = boundary_end_bins(nfull, m)
    sigmain = view(sigmafull, m+1:nfull-m)
    var = boundary_variance(Dspk, G, keep, N, rows, sigmafull[J]) ./ sigmain[rows] .^ 2
    beta = sqrt(maximum(var) / kappa)
    flow = findfirst(keep) - 1
    ell = N / flow
    z = Phiinv(1 - alpha * ell / Ng)
    return beta * sqrt(ell), z, ell, beta
end

function boundary_variance(Dspk::AbstractMatrix, G::AbstractVector, keep::AbstractVector{Bool},
                           N::Integer, rows::AbstractRange, sigmaJ::AbstractVector)
    nJ = size(Dspk, 1); nt = min(Threads.nthreads(), nJ)
    acc = [zeros(length(rows)) for _ in 1:nt]
    Threads.@threads for c in 1:nt
        for a in c:nt:nJ
            h = irfft_cached(whitened_fourier(view(Dspk, a, :), G, keep), N)
            acc[c] .+= view(h, rows) .^ 2 .* sigmaJ[a]^2
        end
    end
    return reduce(+, acc)
end

function tau0_exact(F::SGFilter, Dspk::AbstractMatrix, G::AbstractVector, sigmafull::AbstractVector;
                    alpha::Real = 0.05, lo::Real = 1e-9, hi::Real = 0.5, iters::Integer = 44, drop_nyquist::Bool = true)
    a, b = log(lo), log(hi)
    for _ in 1:iters
        c = (a + b) / 2
        s, z, _, _ = boundary_significance(F, Dspk, G, sigmafull, exp(c); alpha, drop_nyquist)
        s <= z ? (b = c) : (a = c)
    end
    return exp(b)
end

const TAU0_CONSTANTS = Dict(2 => (0.156, 0.419, 1.42), 4 => (0.0501, 0.524, 1.02), 6 => (0.0946, 0.429, 0.791))

boundary_r2(sigmaL::Real, sigmaR::Real) = (sigmaL^2 + sigmaR^2) / (2min(sigmaL, sigmaR)^2)

function tau0_closed(n::Integer, m::Integer, Ng::Integer, r2::Real; alpha::Real = 0.05)
    haskey(TAU0_CONSTANTS, n) || throw(ArgumentError("closed form available for n = 2, 4, 6"))
    A, s, L = TAU0_CONSTANTS[n]
    t = 1e-3
    for _ in 1:3
        ell = L * m * t^(-1 / (2n + 4))
        z = Phiinv(1 - alpha * ell / Ng)
        t = (A^2 * r2 / z^2)^(1 / (2s))
    end
    return t
end
