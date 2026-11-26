function retained_modes(G::AbstractVector{<:Real}, tau::Real; drop_nyquist::Bool = true, N::Integer = 2(length(G) - 1))
    tau >= 0 || throw(ArgumentError("\u03c4 must be non-negative"))
    g2 = abs2.(G); gmax = maximum(g2)
    keep = BitVector(g2 .>= tau * gmax .&& g2 .> 0)
    drop_nyquist && iseven(N) && 2(length(G) - 1) == N && (keep[end] = false)
    return keep
end

n_removed_rfft(keep::AbstractVector{Bool}) = count(!, keep)

removed_modes(keep::AbstractVector{Bool}) = [f - 1 for f in eachindex(keep) if !keep[f]]

mode_weight(f::Integer, N::Integer) = (f == 0 || 2f == N) ? 1 : 2

function retained_fraction(keep::AbstractVector{Bool}, N::Integer)
    s = 0
    for f in eachindex(keep)
        keep[f] && (s += mode_weight(f - 1, N))
    end
    return s / N
end

n_removed(keep::AbstractVector{Bool}, N::Integer) = N - round(Int, retained_fraction(keep, N) * N)

function whitened_fourier(R::AbstractVector{<:Complex}, G::AbstractVector{<:Real}, keep::AbstractVector{Bool})
    C = zeros(ComplexF64, length(R))
    @inbounds for f in eachindex(R)
        keep[f] && (C[f] = R[f] / G[f])
    end
    return C
end

function whiten(r::AbstractVector, G::AbstractVector, keep::AbstractVector{Bool}, sigmahat::AbstractVector)
    N = length(r)
    C = whitened_fourier(rfft(r), G, keep)
    return irfft(C, N) ./ sigmahat
end
