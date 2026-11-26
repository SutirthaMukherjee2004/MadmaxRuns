abstract type LineShape end
struct GaussianLine <: LineShape end
struct LorentzianLine <: LineShape end
struct TopHatLine <: LineShape end
struct AxionLine <: LineShape end

const ALL_SHAPES = (GaussianLine(), LorentzianLine(), TopHatLine(), AxionLine())
shape_name(::GaussianLine) = "Gaussian"
shape_name(::LorentzianLine) = "Lorentzian"
shape_name(::TopHatLine) = "top-hat"
shape_name(::AxionLine) = "axion"
function parse_shape(s::AbstractString)
    t = lowercase(s)
    t in ("gaussian", "gauss") && return GaussianLine()
    t in ("lorentzian", "lorentz") && return LorentzianLine()
    t in ("top-hat", "tophat", "top_hat") && return TopHatLine()
    t in ("axion", "shm", "axion (shm)") && return AxionLine()
    throw(ArgumentError("unknown line shape \"$s\""))
end

fwhm(w::Real) = 2sqrt(2log(2)) * w
Phi_std(z) = 0.5 * erfc(-z / sqrt(2))

const SIGMA_V = 218.0
const V_LAB = 242.1
function axion_speed_cdf(v::Real)
    sv, vl = SIGMA_V, V_LAB
    return sv / (sqrt(2 * pi) * vl) * (exp(-(v + vl)^2 / (2sv^2)) - exp(-(v - vl)^2 / (2sv^2))) +
           0.5 * (erf((v - vl) / (sqrt(2) * sv)) + erf((v + vl) / (sqrt(2) * sv)))
end

struct AxionTable
    u::Vector{Float64}
    c::Vector{Float64}
    upeak::Float64
    fwhm::Float64
end
function build_axion_table(; umax = 6e5, nu = 600_001)
    u = collect(range(0, umax, length = nu)); du = u[2] - u[1]
    c = axion_speed_cdf.(sqrt.(2 .* u))
    pdf = similar(c)
    pdf[2:end-1] = (c[3:end] .- c[1:end-2]) ./ (2du)
    pdf[1] = (c[2] - c[1]) / du; pdf[end] = (c[end] - c[end-1]) / du
    ip = argmax(pdf); h = pdf[ip] / 2
    lo = interp_linear(pdf[1:ip], u[1:ip], h)
    hi = interp_linear(-pdf[ip:end], u[ip:end], -h)
    return AxionTable(u, c, u[ip], hi - lo)
end
function interp_linear(x::AbstractVector, y::AbstractVector, x0::Real)
    x0 <= x[1] && return y[1]
    x0 >= x[end] && return y[end]
    k = searchsortedlast(x, x0)
    return y[k] + (y[k+1] - y[k]) * (x0 - x[k]) / (x[k+1] - x[k])
end
const AXION_TABLE = Ref{Union{Nothing,AxionTable}}(nothing)
const AXION_LOCK = ReentrantLock()
function axion_table()
    t = AXION_TABLE[]
    t === nothing || return t
    lock(AXION_LOCK) do
        AXION_TABLE[] === nothing && (AXION_TABLE[] = build_axion_table())
    end
    return AXION_TABLE[]
end

line_cdf(::GaussianLine, x::Real, w::Real) = Phi_std(x / w)
line_cdf(::LorentzianLine, x::Real, w::Real) = 0.5 + atan(2x / fwhm(w)) / pi
line_cdf(::TopHatLine, x::Real, w::Real) = clamp(x / fwhm(w) + 0.5, 0.0, 1.0)
function line_cdf(::AxionLine, x::Real, w::Real)
    tab = axion_table()
    uu = tab.upeak + x * tab.fwhm / fwhm(w)
    return uu <= 0 ? 0.0 : interp_linear(tab.u, tab.c, uu)
end

function line_profile(shape::LineShape, A::Real, w::Real, nbins::Integer, centre::Real)
    C = [line_cdf(shape, e - centre, w) for e in 0:nbins]
    return A .* diff(C)
end

function line_area_for_peak(shape::LineShape, w::Real, peak::Real; centre_frac::Real = 0.5)
    p = line_profile(shape, 1.0, w, 2001, 1000 + centre_frac)
    return peak / maximum(p)
end
