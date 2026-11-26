struct NullSet
    kind::Symbol
    T::Vector{Float64}
end

function mc_pvalue(T::Real, null::NullSet)
    n = length(null.T)
    k = n - searchsortedfirst(null.T, T) + 1
    return (1 + k) / (1 + n)
end

stream_rng(seed::Integer, i::Integer) = Xoshiro(UInt64(seed) * 0x9e3779b97f4a7c15 + UInt64(i))

function null_iid(S::FilterSetup, n::Integer; seed::Integer = 11)
    T = zeros(n)
    Threads.@threads for i in 1:n
        y = randn(stream_rng(seed, i), S.N)
        Y = rfft(y)
        @inbounds for f in eachindex(Y)
            S.keep[f] || (Y[f] = 0)
        end
        v = irfft(Y, S.N)
        T[i] = squares_T(view(v, S.rows) ./ sqrt(S.kappa))
    end
    return NullSet(:iid, sort!(T))
end

function null_filt(S::FilterSetup, sigmapool::AbstractVector{<:AbstractVector}, n::Integer;
                   seed::Integer = 13, reestimate::Bool = true)
    T = zeros(n); m = S.cfg.m
    Threads.@threads for i in 1:n
        sigma = sigmapool[mod1(i, length(sigmapool))]
        sigmafull = pad_edges(sigma, m)
        x = sigmafull .* randn(stream_rng(seed, i), S.nfull)
        nz = reestimate ? nothing :
             NoiseFit(NaN, NaN, zeros(S.N), collect(sigma), 0, true, false, Float64[], Float64[], :given, S.cfg.W)
        T[i] = run_chain(S, x; noise = nz).T
    end
    return NullSet(:filt, sort!(T))
end
