Base.@kwdef mutable struct StudyOptions
    n::Int = 4
    m::Int = 64
    kind::Symbol = :SG
    nspec::Int = 5000
    nnull::Int = 5000
    nfilt::Int = 5000
    seed::Int = 2026
    outdir::String = "results"
    levels::Vector{Float64} = [0.0012, 0.01, 0.05, 0.173]
    models::Vector{NoiseModel} = [Homoscedastic, Heteroscedastic]
    tau::Float64 = 1e-7
    W::Int = 96
    filt_reestimate::Bool = false
    taus::Vector{Float64} = Float64[]
    tau_min::Float64 = NaN
    tau0_homo::Float64 = NaN
    tau0_hetero::Float64 = NaN
    shapes::Vector{String} = ["Gaussian", "Lorentzian", "top-hat", "axion"]
    widths::Vector{Float64} = collect(2.0:2.0:12.0)
    amplitudes::Vector{Float64} = collect(0.0:4.0:28.0)
    level::Float64 = 0.05
    quick::Bool = false
end

base_config(o::StudyOptions; kw...) = with(PipelineConfig(n = o.n, m = o.m, kind = o.kind, tau = o.tau, W = o.W); kw...)
level_tag(lv) = replace(@sprintf("%.2f", 100lv), "." => "p") * "pct"
result_path(o::StudyOptions, name) = joinpath(o.outdir, name * ".jls")

paper_tau_grid() = sort(unique(vcat(exp10.(range(-7, log10(0.5), length = 32)), exp10.(range(log10(3e-6), log10(3e-3), length = 31)))))
tau_grid(o::StudyOptions) = isempty(o.taus) ? paper_tau_grid() : sort(o.taus)

spectrum_rng(o::StudyOptions, case::Integer, s::Integer) = stream_rng(o.seed + 1_000_003 * case, s)

function minimum_threshold(o::StudyOptions)
    tmin = o.tau_min
    if isnan(tmin) && isfile(result_path(o, "tauscan"))
        tmin = load_result(result_path(o, "tauscan")).tau_min
    end
    return isnan(tmin) ? 5.4e-5 : tmin
end

function thresholds(o::StudyOptions)
    tmin = minimum_threshold(o); t0h = o.tau0_homo; t0x = o.tau0_hetero
    if (isnan(t0h) || isnan(t0x)) && isfile(result_path(o, "tau0"))
        r = load_result(result_path(o, "tau0"))
        isnan(t0h) && (t0h = r.tau0[(o.n, o.m, "homo")])
        isnan(t0x) && (t0x = r.tau0[(o.n, o.m, "hetero")])
    end
    isnan(t0h) && (t0h = 2.9e-4)
    isnan(t0x) && (t0x = 3.9e-4)
    return tmin, t0h, t0x
end

struct SpectrumState
    noise::NoiseFit
    r::Vector{Float64}
    R::Vector{ComplexF64}
end
function spectrum_state(S::FilterSetup, x::AbstractVector; noise = nothing)
    nf = noise === nothing ? estimate_noise(S.F, S.wsg, x; W = S.cfg.W, estimator = S.cfg.estimator,
                                            tol = S.cfg.tol, maxit = S.cfg.maxit) : noise
    validate_noise_fit(nf, S.N)
    r = linear_residual(S.F, x)
    return SpectrumState(nf, r, rfft(r))
end

function tau_whitened(S::FilterSetup, sigma::AbstractVector, R::AbstractVector)
    length(sigma) == S.N || throw(DimensionMismatch("expected $(S.N) noise bins"))
    all(s -> isfinite(s) && s > 0, sigma) || error("whitening requires finite positive noise scales")
    C = whitened_fourier(R, S.G, S.keep)
    c = irfft(C, S.N)
    z = c[S.rows] ./ sigma[S.rows] ./ sqrt(S.kappa)
    all(isfinite, z) || error("whitened residual contains non-finite values")
    return z
end

function tau_step(S::FilterSetup, sigma::AbstractVector, R::AbstractVector; removal::Bool = true)
    z = tau_whitened(S, sigma, R)
    T = squares_T(z)
    isfinite(T) || error("non-finite SQUARES statistic")
    Tp = NaN; zp = z
    if removal && size(S.H, 1) > 0
        B = boundary_patterns(S.H, pad_edges(sigma, S.cfg.m)[S.J], sigma[S.rows], S.kappa; K = S.cfg.nboundary)
        zp = remove_boundary(z, B)
        all(isfinite, zp) || error("boundary-projected residual contains non-finite values")
        Tp = squares_T(zp)
        isfinite(Tp) || error("non-finite boundary-projected SQUARES statistic")
    end
    return (; z, zp, T, Tp)
end
tau_step(S::FilterSetup, st::SpectrumState; removal::Bool = true) = tau_step(S, st.noise.sigma, st.R; removal)

struct SpectrumSet
    name::String
    gen::Function
    sigma::Vector{Vector{Float64}}
end
Base.length(s::SpectrumSet) = length(s.sigma)

function fit_set(S::FilterSetup, name::AbstractString, gen::Function, n::Integer)
    sigma = Vector{Vector{Float64}}(undef, n)
    Threads.@threads for i in 1:n
        cfg = S.cfg
        nf = estimate_noise(S.F, S.wsg, gen(i); W = cfg.W, estimator = cfg.estimator, tol = cfg.tol, maxit = cfg.maxit)
        try
            validate_noise_fit(nf, S.N)
        catch err
            throw(ErrorException("$name spectrum $i: $(sprint(showerror, err))"))
        end
        sigma[i] = nf.sigma
    end
    return SpectrumSet(String(name), gen, sigma)
end

function data_generator(o::StudyOptions, cfg::PipelineConfig, case::Integer; model::NoiseModel, level::Real,
                        signal = nothing, background_on::Bool = true)
    bmed = median_background(; grid = cfg.grid)
    return i -> simulate_spectrum(spectrum_rng(o, case, i), cfg; model, level, bmed, signal, background_on)[1]
end

function filt_generator(sigmapool::AbstractVector, m::Integer, seed::Integer)
    return i -> pad_edges(Float64.(sigmapool[mod1(i, length(sigmapool))]), m) .* randn(stream_rng(seed, i), length(sigmapool[1]) + 2m)
end

case_id(model::NoiseModel, level::Real; background_on::Bool = true) =
    100_000 * (Int(model) + 1) + round(Int, 1e4 * level) + (background_on ? 0 : 50_000)

function eval_set(S::FilterSetup, set::SpectrumSet; removal::Bool = true, given_sigma = nothing)
    n = length(set); T = zeros(n); Tp = fill(NaN, n)
    Threads.@threads for i in 1:n
        x = set.gen(i)
        R = rfft(linear_residual(S.F, x))
        sigma = given_sigma === nothing ? Float64.(set.sigma[i]) : Float64.(given_sigma[mod1(i, length(given_sigma))])
        st = tau_step(S, sigma, R; removal)
        T[i] = st.T; Tp[i] = st.Tp
    end
    return T, Tp
end

function whitened_line(S::FilterSetup, s::AbstractVector, sigma::AbstractVector)
    C = whitened_fourier(rfft(linear_residual(S.F, s)), S.G, S.keep)
    return (irfft(C, S.N) ./ sigma)[S.rows] ./ sqrt(S.kappa)
end

axion_offset(::SpectrumGrid = PAPER_GRID) = 10.5e6

function axion_signal(A::Union{Nothing,Real} = nothing; w::Real = 3.4)
    shape = AxionLine()
    area = A === nothing ? line_area_for_peak(shape, w, 3.0) : Float64(A)
    return shape, area, w, axion_offset()
end

function tau_setups(cfg::PipelineConfig, taus; with_pk = true, with_boundary = true)
    F = SGFilter(cfg.n, cfg.m; kind = cfg.kind)
    D = with_boundary ? spike_fourier(F, cfg.grid.nbins) : nothing
    return [filter_setup(with_tau(cfg, tau); with_pk, with_boundary, Dspk = D) for tau in taus]
end

fraction_below(p, alpha) = count(<=(alpha), p) / length(p)

function p_histogram(p::AbstractVector; nb::Integer = 20)
    h = zeros(Int, nb)
    for v in p
        isfinite(v) && 0 <= v <= 1 || throw(ArgumentError("p-value histogram received an invalid probability"))
        h[clamp(floor(Int, v * nb) + 1, 1, nb)] += 1
    end
    return h
end
