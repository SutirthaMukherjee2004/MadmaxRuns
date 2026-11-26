struct FilterSetup
    cfg::PipelineConfig
    F::SGFilter
    wsg::WSGCache
    nfull::Int
    N::Int
    rows::UnitRange{Int}
    Ng::Int
    G::Vector{Float64}
    keep::BitVector
    K::Vector{Int}
    kappa::Float64
    q::Float64
    Dspk::Matrix{ComplexF64}
    H::Matrix{Float64}
    J::Vector{Int}
    pk::PKModel
    pktab::Union{Nothing,PKTable}
    pbctab::Union{Nothing,PKTable}
end

function filter_setup(cfg::PipelineConfig; with_pk::Bool = true, with_boundary::Bool = cfg.boundary_removal,
                      Dspk::Union{Nothing,Matrix{ComplexF64}} = nothing, pk_tmax::Real = 150.0)
    cfg.grid.nbins > 4cfg.m || throw(ArgumentError("the spectrum must have more than 4m bins"))
    isfinite(cfg.grid.df_hz) && cfg.grid.df_hz > 0 && isfinite(cfg.grid.fref_hz) ||
        throw(ArgumentError("the frequency grid must have a finite origin and positive finite spacing"))
    isfinite(cfg.tau) && 0 <= cfg.tau <= 1 || throw(ArgumentError("tau must be in [0, 1]"))
    cfg.nboundary >= 0 || throw(ArgumentError("nboundary must be nonnegative"))
    F = SGFilter(cfg.n, cfg.m; kind = cfg.kind)
    nfull = cfg.grid.nbins; N = nfull - 2cfg.m
    rows = tested_bins(N, cfg.m); Ng = length(rows)
    G = transfer_function(F, N)
    all(isfinite, G) || error("filter transfer function contains non-finite values")
    keep = retained_modes(G, cfg.tau; drop_nyquist = cfg.drop_nyquist, N)
    K = removed_modes(keep)
    kappa = retained_fraction(keep, N)
    kappa > 0 || throw(ArgumentError("the Fourier mask retains no modes"))
    Dspk = with_boundary ? (Dspk === nothing ? spike_fourier(F, nfull) : Dspk) : zeros(ComplexF64, 0, 0)
    if with_boundary
        size(Dspk) == (4cfg.m, length(G)) || throw(DimensionMismatch("boundary spike response dimensions are invalid"))
        all(isfinite, Dspk) || error("boundary spike response contains non-finite values")
    end
    H = with_boundary ? boundary_responses(Dspk, G, keep, N, rows) : zeros(0, 0)
    all(isfinite, H) || error("whitened boundary response contains non-finite values")
    pk = PKModel(N, Ng, K; qmax = cfg.qmax, order = cfg.pk_order)
    pktab = with_pk ? pk_table(pk; tmax = pk_tmax) : nothing
    pbctab = with_pk ? pbc_table(Ng; tmax = pk_tmax) : nothing
    return FilterSetup(cfg, F, WSGCache(F, nfull), nfull, N, rows, Ng, G, keep, K, kappa,
                       noise_pass_fraction(F), Dspk, H, boundary_end_bins(nfull, cfg.m), pk, pktab, pbctab)
end
