# Thread-safe real FFTs with plans that are created once per length and then reused.
#
# FFTW.rfft / FFTW.irfft create a new plan on every call, and plan creation holds a global lock in
# FFTW. When many spectra are processed in parallel this lock serialises the threads. Executing an
# existing plan on new arrays is thread-safe in FFTW, so the plans are cached here by length. They are
# made with FFTW.ESTIMATE (as rfft/irfft do) and FFTW.UNALIGNED, so any Vector can be passed. The
# transforms are the same discrete Fourier transforms; results agree with rfft/irfft to rounding.
# The plans are not stored in FilterSetup, which is serialised into checkpoints.

const FFT_PLAN_LOCK = ReentrantLock()
const RFFT_PLANS = Dict{Int,Any}()
const BRFFT_PLANS = Dict{Int,Any}()

function cached_plan(cache::Dict{Int,Any}, n::Integer, make)
    p = lock(FFT_PLAN_LOCK) do
        get(cache, n, nothing)
    end
    p === nothing || return p
    return lock(FFT_PLAN_LOCK) do
        get!(cache, n) do
            make()
        end
    end
end

function rfft_cached(x::AbstractVector{<:Real})
    n = length(x)
    v = x isa Vector{Float64} ? x : Vector{Float64}(x)
    p = cached_plan(RFFT_PLANS, n, () -> plan_rfft(zeros(n); flags = FFTW.ESTIMATE | FFTW.UNALIGNED))
    y = Vector{ComplexF64}(undef, div(n, 2) + 1)
    mul!(y, p, v)
    return y
end

function irfft_cached(X::AbstractVector{<:Complex}, n::Integer)
    length(X) == div(n, 2) + 1 || throw(DimensionMismatch("irfft of length $n needs $(div(n, 2) + 1) coefficients"))
    p = cached_plan(BRFFT_PLANS, n, () -> plan_brfft(zeros(ComplexF64, div(n, 2) + 1), n; flags = FFTW.ESTIMATE | FFTW.UNALIGNED))
    y = Vector{Float64}(undef, n)
    mul!(y, p, Vector{ComplexF64}(X))       # c2r transforms overwrite their input: always pass a copy
    rmul!(y, 1 / n)
    return y
end

# In-place inverse real FFT into y; X is used as work space and is overwritten.
function irfft_cached!(y::Vector{Float64}, X::Vector{ComplexF64}, n::Integer)
    length(y) == n && length(X) == div(n, 2) + 1 || throw(DimensionMismatch("irfft! buffers do not match length $n"))
    p = cached_plan(BRFFT_PLANS, n, () -> plan_brfft(zeros(ComplexF64, div(n, 2) + 1), n; flags = FFTW.ESTIMATE | FFTW.UNALIGNED))
    mul!(y, p, X)
    rmul!(y, 1 / n)
    return y
end

# Solves A x = y in place for a symmetric positive definite A (full matrix given) by Cholesky
# factorisation, written out so that the many small solves of the weighted SG filter and the window
# fit do not go through LAPACK, whose small calls serialise across threads. Returns y.
function spd_solve!(A::AbstractMatrix{Float64}, y::AbstractVector{Float64})
    n = size(A, 1)
    @inbounds for j in 1:n
        s = A[j, j]
        for k in 1:j-1
            s -= A[j, k]^2
        end
        s > 0 || throw(PosDefException(j))
        d = sqrt(s); A[j, j] = d
        for i in j+1:n
            v = A[i, j]
            for k in 1:j-1
                v -= A[i, k] * A[j, k]
            end
            A[i, j] = v / d
        end
    end
    @inbounds for i in 1:n
        v = y[i]
        for k in 1:i-1
            v -= A[i, k] * y[k]
        end
        y[i] = v / A[i, i]
    end
    @inbounds for i in n:-1:1
        v = y[i]
        for k in i+1:n
            v -= A[k, i] * y[k]
        end
        y[i] = v / A[i, i]
    end
    return y
end
