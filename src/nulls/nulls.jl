struct NullSet
    kind::Symbol
    T::Vector{Float64}
end

function mc_pvalue(T::Real, null::NullSet)
    n = length(null.T)
    n > 0 || throw(ArgumentError("Monte Carlo p-value requires at least one null replica"))
    isfinite(T) || throw(ArgumentError("Monte Carlo p-value requires a finite SQUARES statistic"))
    k = n - searchsortedfirst(null.T, T) + 1
    return (1 + k) / (1 + n)
end

stream_rng(seed::Integer, i::Integer) = Xoshiro(UInt64(seed) * 0x9e3779b97f4a7c15 + UInt64(i))

function null_iid(S::FilterSetup, n::Integer; seed::Integer = 11, first_replica::Integer = 1)
    n > 0 || throw(ArgumentError("n must be positive"))
    first_replica > 0 || throw(ArgumentError("first_replica must be positive"))
    T = zeros(n)
    Threads.@threads for i in 1:n
        replica = first_replica + i - 1
        y = randn(stream_rng(seed, replica), S.N)
        Y = rfft_cached(y)
        @inbounds for f in eachindex(Y)
            S.keep[f] || (Y[f] = 0)
        end
        v = irfft_cached(Y, S.N)
        T[i] = first(squares_T_run(view(v, S.rows) ./ sqrt(S.kappa)))
        isfinite(T[i]) || error("non-finite SQUARES statistic in i.i.d. null replica $replica")
    end
    return NullSet(:iid, sort!(T))
end

function null_filt(S::FilterSetup, sigmapool::AbstractVector{<:AbstractVector}, n::Integer;
                   seed::Integer = 13, reestimate::Bool = false, first_replica::Integer = 1)
    n > 0 || throw(ArgumentError("n must be positive"))
    first_replica > 0 || throw(ArgumentError("first_replica must be positive"))
    isempty(sigmapool) && throw(ArgumentError("filtered null requires at least one fitted noise profile"))
    for sigma in sigmapool
        length(sigma) == S.N || throw(DimensionMismatch("expected $(S.N) fitted noise bins"))
        all(s -> isfinite(s) && s > 0, sigma) || throw(ArgumentError("noise profile must be finite and positive"))
    end
    T = zeros(n); m = S.cfg.m
    Threads.@threads for i in 1:n
        replica = first_replica + i - 1
        sigma = sigmapool[mod1(replica, length(sigmapool))]
        sigmafull = pad_edges(sigma, m)
        x = sigmafull .* randn(stream_rng(seed, replica), S.nfull)
        nz = reestimate ? nothing :
             NoiseFit(NaN, NaN, zeros(S.N), collect(sigma), 0, true, false, Float64[], Float64[], :given, S.cfg.W)
        T[i] = run_chain(S, x; noise = nz, compute_pvalues = false).T
    end
    return NullSet(:filt, sort!(T))
end
