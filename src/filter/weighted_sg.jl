struct WSGCache
    F::SGFilter
    nfull::Int
    ckern::Vector{Vector{ComplexF64}}
end
function WSGCache(F::SGFilter, nfull::Integer)
    m, n = F.m, F.n
    fw = F.kind == :SGW ? [cos(pi * d / (2(m + 1)))^4 for d in -m:m] : ones(2m + 1)
    d = collect(-m:m) ./ m
    ck = Vector{Vector{ComplexF64}}(undef, 2n + 1)
    for p in 0:2n
        k = zeros(nfull); k[1:2m+1] = fw .* d .^ p
        ck[p+1] = conj.(rfft_cached(k))
    end
    return WSGCache(F, nfull, ck)
end

function _corr(Y::Vector{ComplexF64}, ck::Vector{ComplexF64}, nfull::Int, N::Int)
    v = irfft_cached(Y .* ck, nfull)
    return v[1:N]
end

function weighted_smooth(c::WSGCache, x::AbstractVector, wt::AbstractVector)
    F = c.F; m, n, nf = F.m, F.n, c.nfull
    length(x) == nf == length(wt) || throw(DimensionMismatch("x and wt need $(nf) bins"))
    N = nf - 2m
    sw = nf / sum(wt)
    w = wt .* sw
    W = rfft_cached(w); WX = rfft_cached(w .* x)
    # Sliding weighted moments sum_d w_(i+d) d^p (p = 0..2n) and sum_d w_(i+d) x_(i+d) d^p (p = 0..n),
    # as FFT correlations with the kernels d^p; columns of Mom and Yv. Buffers are reused.
    Mom = Matrix{Float64}(undef, N, 2n + 1); Yv = Matrix{Float64}(undef, N, n + 1)
    Ybuf = similar(W); vbuf = Vector{Float64}(undef, nf)
    for p in 0:2n
        Ybuf .= W .* c.ckern[p+1]
        irfft_cached!(vbuf, Ybuf, nf)
        copyto!(view(Mom, :, p + 1), view(vbuf, 1:N))
    end
    for p in 0:n
        Ybuf .= WX .* c.ckern[p+1]
        irfft_cached!(vbuf, Ybuf, nf)
        copyto!(view(Yv, :, p + 1), view(vbuf, 1:N))
    end
    P = n + 1
    out = zeros(N)
    A = zeros(P, P); y = zeros(P)
    @inbounds for i in 1:N
        for a in 1:P, b in 1:P
            A[a, b] = Mom[i, a+b-1]
        end
        for a in 1:P
            y[a] = Yv[i, a]
        end
        spd_solve!(A, y)
        out[i] = y[1]
    end
    return out
end

function weighted_smooth_direct(F::SGFilter, x::AbstractVector, wt::AbstractVector)
    m = F.m; N = length(x) - 2m
    V = sg_vandermonde(F.n, m)
    fw = F.kind == :SGW ? [cos(pi * d / (2(m + 1)))^4 for d in -m:m] : ones(2m + 1)
    return [begin
        j = i:i+2m
        Wd = Diagonal(fw .* wt[j])
        ((V' * Wd * V) \ (V' * Wd * x[j]))[1]
    end for i in 1:N]
end
