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
        ck[p+1] = conj.(rfft(k))
    end
    return WSGCache(F, nfull, ck)
end

function _corr(Y::Vector{ComplexF64}, ck::Vector{ComplexF64}, nfull::Int, N::Int)
    v = irfft(Y .* ck, nfull)
    return v[1:N]
end

function weighted_smooth(c::WSGCache, x::AbstractVector, wt::AbstractVector)
    F = c.F; m, n, nf = F.m, F.n, c.nfull
    length(x) == nf == length(wt) || throw(DimensionMismatch("x and wt need $(nf) bins"))
    N = nf - 2m
    sw = nf / sum(wt)
    w = wt .* sw
    W = rfft(w); WX = rfft(w .* x)
    Mom = [_corr(W, c.ckern[p+1], nf, N) for p in 0:2n]
    Yv = [_corr(WX, c.ckern[p+1], nf, N) for p in 0:n]
    P = n + 1
    out = zeros(N)
    A = zeros(P, P); y = zeros(P)
    @inbounds for i in 1:N
        for a in 1:P, b in 1:P
            A[a, b] = Mom[a+b-1][i]
        end
        for a in 1:P
            y[a] = Yv[a][i]
        end
        ldiv!(cholesky!(Symmetric(A)), y)
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
