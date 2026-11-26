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

function validate_noise_fit(fit::NoiseFit, n::Integer)
    fit.converged || throw(ArgumentError("noise fit did not converge after $(fit.iterations) iterations"))
    length(fit.sigma) == length(fit.bhat) == n ||
        throw(DimensionMismatch("noise fit and background must have $n bins"))
    all(v -> isfinite(v) && v > 0, fit.sigma) ||
        throw(ArgumentError("fitted noise must be finite and strictly positive in every bin"))
    all(isfinite, fit.bhat) || throw(ArgumentError("fitted background contains a non-finite value"))
    if fit.estimator != :given
        all(isfinite, (fit.sigma0, fit.sigma1)) || throw(ArgumentError("noise fit has non-finite parameters"))
    end
    return fit
end

function chain_noise(S, x, noise, checkpoint)
    cfg = S.cfg
    nf = noise === nothing ?
         estimate_noise(S.F, S.wsg, x; W = cfg.W, estimator = cfg.estimator, tol = cfg.tol, maxit = cfg.maxit) :
         noise
    try
        return validate_noise_fit(nf, S.N)
    catch
        if checkpoint !== nothing
            try
                path = write_checkpoint(checkpoint, :noise_fit, nf; failed = true)
                log_message("Failed noise-fit diagnostics saved to $path")
            catch err
                log_message("Could not save failed noise-fit diagnostics: $(sprint(showerror, err))")
            end
        end
        rethrow()
    end
end

function chain_whitening(S, x, sigma)
    r = linear_residual(S.F, x)
    rho = r ./ sigma
    R = rfft_cached(r)
    C = whitened_fourier(R, S.G, S.keep)
    c = irfft_cached(C, S.N) ./ sigma
    z = c[S.rows] ./ sqrt(S.kappa)
    all(v -> all(isfinite, v), (r, rho, R, C, c, z)) || error("whitening produced non-finite values")
    return (; r, rho, R, C, c, z)
end

function chain_projection(S, sigma, z)
    S.cfg.boundary_removal || return (; zp = z, boundary = nothing)
    sigmafull = pad_edges(sigma, S.cfg.m)
    B = boundary_patterns(S.H, sigmafull[S.J], sigma[S.rows], S.kappa; K = S.cfg.nboundary)
    all(v -> all(isfinite, v), (B.s, B.V, B.uncaptured)) ||
        error("boundary basis contains non-finite values")
    zp = remove_boundary(z, B)
    all(isfinite, zp) || error("boundary projection produced non-finite values")
    return (; zp, boundary = B)
end

# T is the largest chi^2 of any success run, Eq. (2). It is found by one linear scan
# (squares_T_run), which gives the same value as RunStatistics.t_obs bit for bit. For the analysed
# spectrum itself (verify = true) both are computed and must agree, and the RunStatistics value is
# reported; Monte Carlo replicas use the scan alone.
function chain_statistics(S, residual, sigma, zp; verify::Bool = false)
    T, i0, i1 = squares_T_run(zp)
    zr = residual[S.rows] ./ (sigma[S.rows] .* sqrt(S.q))
    all(isfinite, zr) || error("normalized raw residual contains non-finite values")
    Tr = first(squares_T_run(zr))
    if verify
        Tref = squares_T(zp); Trref = squares_T(zr)
        isapprox(T, Tref; rtol = 1e-12) && isapprox(Tr, Trref; rtol = 1e-12) ||
            error("SQUARES T from the run scan ($T, $Tr) disagrees with RunStatistics.t_obs ($Tref, $Trref)")
        T, Tr = Tref, Trref
    end
    all(v -> isfinite(v) && v >= 0, (T, Tr)) || error("SQUARES statistic is invalid")
    return (; T, run = (i0, i1), Tr)
end

function chain_pvalues(S, T, Tr)
    p = tabulated_pvalue(t -> pk_pvalue(S.pk, t), S.pktab, T)
    pr = tabulated_pvalue(t -> pbc_pvalue(t, S.Ng), S.pbctab, Tr)
    all(v -> isfinite(v) && 0 < v <= 1, (p, pr)) || error("p-value calculation produced an invalid probability")
    return (; p, pr)
end

function run_chain(S::FilterSetup, x::AbstractVector; noise::Union{Nothing,NoiseFit} = nothing,
                   compute_pvalues::Bool = true)
    return _run_chain(S, x, nothing; noise, compute_pvalues)
end

function _run_chain(S::FilterSetup, x::AbstractVector, checkpoint;
                    noise::Union{Nothing,NoiseFit} = nothing, compute_pvalues::Bool = true)
    length(x) == S.nfull || throw(DimensionMismatch("spectrum must have $(S.nfull) bins"))
    all(isfinite, x) || throw(ArgumentError("spectrum contains non-finite values"))
    nf = run_stage(checkpoint, :noise_fit) do
        chain_noise(S, x, noise, checkpoint)
    end
    sigma = nf.sigma
    w = run_stage(checkpoint, :whitening) do
        chain_whitening(S, x, sigma)
    end
    b = run_stage(checkpoint, :boundary_projection) do
        chain_projection(S, sigma, w.z)
    end
    obs = run_stage(checkpoint, :statistics) do
        chain_statistics(S, w.r, sigma, b.zp; verify = compute_pvalues)
    end
    values = if compute_pvalues
        run_stage(checkpoint, :pvalues) do
            chain_pvalues(S, obs.T, obs.Tr)
        end
    else
        (; p = NaN, pr = NaN)
    end
    return ChainResult(nf, w.r, w.rho, w.R, w.C, w.c, w.z, b.zp, b.boundary,
                       obs.T, obs.run, values.p, obs.Tr, values.pr)
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
