struct SGFilter
    n::Int
    m::Int
    kind::Symbol
    S::Vector{Float64}
    g::Vector{Float64}
end
function SGFilter(n::Integer, m::Integer; kind::Symbol = :SG)
    n >= 0 && m >= 1 && 2m + 1 > n || throw(ArgumentError("need 2m+1 > n \u2265 0"))
    kind in (:SG, :SGW) || throw(ArgumentError("kind must be :SG or :SGW"))
    S = kind == :SG ? sg_kernel(n, m) : sgw_kernel(n, m)
    g = -copy(S); g[m+1] += 1
    return SGFilter(n, m, kind, S, g)
end
Base.show(io::IO, F::SGFilter) = print(io, "$(F.kind)(n=$(F.n), m=$(F.m))")

sg_vandermonde(n, m) = hcat([(collect(-m:m) ./ m) .^ p for p in 0:n]...)

function sg_kernel(n::Integer, m::Integer)
    V = sg_vandermonde(n, m)
    e = zeros(n + 1); e[1] = 1
    return V * ((V' * V) \ e)
end

function sgw_kernel(n::Integer, m::Integer)
    V = sg_vandermonde(n, m)
    w = [cos(pi * d / (2(m + 1)))^4 for d in -m:m]
    e = zeros(n + 1); e[1] = 1
    return w .* (V * ((V' * Diagonal(w) * V) \ e))
end

residual_kernel(F::SGFilter) = F.g

noise_pass_fraction(F::SGFilter) = sum(abs2, F.g)

kernel_moments(F::SGFilter, p::Integer) = sum(F.g[k] * (-(k - F.m - 1))^p for k in eachindex(F.g))

function smoothing_spread(F::SGFilter)
    q = noise_pass_fraction(F)
    return sum(F.S[k]^2 * (k - F.m - 1)^2 for k in eachindex(F.S) if k != F.m + 1) / q
end

function transfer_function(F::SGFilter, N::Integer)
    N > 2F.m || throw(ArgumentError("N must exceed 2m"))
    k = zeros(N)
    for (j, d) in enumerate(-F.m:F.m)
        k[mod(d, N)+1] += F.g[j]
    end
    return real.(rfft(k))
end

function linear_residual!(r::AbstractVector, F::SGFilter, x::AbstractVector)
    m = F.m; N = length(x) - 2m
    length(r) == N || throw(DimensionMismatch("r must have length(x) - 2m"))
    S = F.S
    @inbounds for i in 1:N
        s = 0.0
        @simd for k in 1:2m+1
            s += S[k] * x[i+k-1]
        end
        r[i] = x[i+m] - s
    end
    return r
end
linear_residual(F::SGFilter, x::AbstractVector) = linear_residual!(zeros(length(x) - 2F.m), F, x)

function smooth_unweighted(F::SGFilter, x::AbstractVector)
    m = F.m
    return x[m+1:end-m] .- linear_residual(F, x)
end
