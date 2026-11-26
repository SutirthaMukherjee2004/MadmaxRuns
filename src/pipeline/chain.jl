struct ChainResult
    noise::NoiseFit
    r::Vector{Float64}
    rho::Vector{Float64}
    R::Vector{ComplexF64}
    C::Vector{ComplexF64}
    c::Vector{Float64}
    z::Vector{Float64}
    zp::Vector{Float64}
    boundary::Union{Nothing,BoundaryBasis}
    T::Float64
    run::Tuple{Int,Int}
    p_pk::Float64
    T_raw::Float64
    p_raw::Float64
end

function run_chain(S::FilterSetup, x::AbstractVector; noise::Union{Nothing,NoiseFit} = nothing)
    cfg = S.cfg; m = cfg.m
    length(x) == S.nfull || throw(DimensionMismatch("spectrum must have $(S.nfull) bins"))
    nf = noise === nothing ?
         estimate_noise(S.F, S.wsg, x; W = cfg.W, estimator = cfg.estimator, tol = cfg.tol, maxit = cfg.maxit) :
         noise
    sigma = nf.sigma
    r = linear_residual(S.F, x)
    rho = r ./ sigma
    R = rfft(r)
    C = whitened_fourier(R, S.G, S.keep)
    c = irfft(C, S.N) ./ sigma
    z = c[S.rows] ./ sqrt(S.kappa)
    B = nothing; zp = z
    if cfg.boundary_removal
        sigmafull = pad_edges(sigma, m)
        B = boundary_patterns(S.H, sigmafull[S.J], sigma[S.rows], S.kappa; K = cfg.nboundary)
        zp = remove_boundary(z, B)
    end
    T = squares_T(zp)
    _, i0, i1 = squares_T_run(zp)
    p = pvalue(S.pktab, T)
    zr = r[S.rows] ./ (sigma[S.rows] .* sqrt(S.q))
    Tr = squares_T(zr); pr = pvalue(S.pbctab, Tr)
    return ChainResult(nf, r, rho, R, C, c, z, zp, B, T, (i0, i1), p, Tr, pr)
end

function simulate_spectrum(rng::AbstractRNG, cfg::PipelineConfig; model::NoiseModel, level::Real,
                           bmed::Real = median_background(; grid = cfg.grid), signal = nothing,
                           background_on::Bool = true)
    g = cfg.grid
    b = background_on ? background(rng; grid = g) : zeros(g.nbins)
    sigma0, sigma1 = noise_parameters(model, level, bmed)
    bsigma = background_on ? b : background(rng; grid = g)
    sigma = noise_sigma(sigma0, sigma1, bsigma)
    x = b .+ sigma .* randn(rng, g.nbins)
    s = zeros(g.nbins)
    if signal !== nothing
        shape, A, w, off = signal
        s = sigma .* line_profile(shape, A, w, g.nbins, frequency_bin(g, off) - 0.5)
        x .+= s
    end
    return x, b, sigma, s
end

function inject_line!(x::AbstractVector, sigma::AbstractVector, shape::LineShape, A::Real, w::Real, centre_bin::Real)
    x .+= sigma .* line_profile(shape, A, w, length(x), centre_bin - 0.5)
    return x
end

function analyse_spectrum(S::FilterSetup, x::AbstractVector; iid = nothing, filt = nothing)
    res = run_chain(S, x)
    p_iid = iid === nothing ? NaN : mc_pvalue(res.T, iid)
    p_filt = filt === nothing ? NaN : mc_pvalue(res.T, filt)
    return (; result = res, T = res.T, p_pk = res.p_pk, p_iid_mc = p_iid, p_filt_mc = p_filt)
end
