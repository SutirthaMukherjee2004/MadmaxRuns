boundary_end_bins(nfull::Integer, m::Integer) = vcat(1:2m, nfull-2m+1:nfull)

tested_bins(N::Integer, m::Integer) = m+1:N-m

function spike_fourier(F::SGFilter, nfull::Integer)
    J = boundary_end_bins(nfull, F.m); N = nfull - 2F.m
    D = zeros(ComplexF64, length(J), div(N, 2) + 1)
    Threads.@threads for a in eachindex(J)
        x = zeros(nfull); x[J[a]] = 1.0
        D[a, :] = rfft(linear_residual(F, x))
    end
    return D
end

function boundary_responses(Dspk::AbstractMatrix, G::AbstractVector, keep::AbstractVector{Bool},
                            N::Integer, rows::AbstractRange = 1:N)
    nJ = size(Dspk, 1)
    H = zeros(length(rows), nJ)
    Threads.@threads for a in 1:nJ
        C = whitened_fourier(view(Dspk, a, :), G, keep)
        h = irfft(C, N)
        H[:, a] = view(h, rows)
    end
    return H
end

struct BoundaryBasis
    s::Vector{Float64}
    V::Matrix{Float64}
    uncaptured::Vector{Float64}
end

function boundary_patterns(H::AbstractMatrix, sigmaJ::AbstractVector, sigmatest::AbstractVector, kappa::Real; K::Integer = 6)
    nJ = size(H, 2); a = sigmaJ ./ sqrt(kappa)
    Gm = zeros(nJ, nJ); blk = 2048
    buf = zeros(min(blk, size(H, 1)), nJ)
    for r0 in 1:blk:size(H, 1)
        r = r0:min(r0 + blk - 1, size(H, 1)); B = view(buf, 1:length(r), :)
        B .= view(H, r, :) ./ view(sigmatest, r)
        BLAS.syrk!('U', 'T', 1.0, B, 1.0, Gm)
    end
    Gm .*= a .* a'
    E = eigen(Symmetric(Gm, :U))
    ord = sortperm(E.values; rev = true)
    lambda = max.(E.values[ord], 0.0); U = E.vectors[:, ord]
    s = sqrt.(lambda)
    K = min(K, count(>(0), lambda))
    V = (H * (a .* U[:, 1:K])) ./ sigmatest
    for k in 1:K
        V[:, k] ./= s[k]
    end
    tot = sum(abs2, s)
    unc = iszero(tot) ? zeros(length(s)) :
          [max(0.0, 1 - sum(abs2, s[1:k]) / tot) for k in eachindex(s)]
    return BoundaryBasis(s, V, unc)
end

remove_boundary(z::AbstractVector, B::BoundaryBasis) = z .- B.V * (B.V' * z)
