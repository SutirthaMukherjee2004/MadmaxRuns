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
    pktab::PKTable
    pbctab::PKTable
end

function filter_setup(cfg::PipelineConfig; with_pk::Bool = true, with_boundary::Bool = cfg.boundary_removal,
                      Dspk::Union{Nothing,Matrix{ComplexF64}} = nothing, pk_tmax::Real = 150.0)
    F = SGFilter(cfg.n, cfg.m; kind = cfg.kind)
    nfull = cfg.grid.nbins; N = nfull - 2cfg.m
    rows = tested_bins(N, cfg.m); Ng = length(rows)
    G = transfer_function(F, N)
    keep = retained_modes(G, cfg.tau; drop_nyquist = cfg.drop_nyquist, N)
    K = removed_modes(keep)
    kappa = retained_fraction(keep, N)
    Dspk = with_boundary ? (Dspk === nothing ? spike_fourier(F, nfull) : Dspk) : zeros(ComplexF64, 0, 0)
    H = with_boundary ? boundary_responses(Dspk, G, keep, N, rows) : zeros(0, 0)
    pk = PKModel(N, Ng, K; qmax = cfg.qmax, order = cfg.pk_order)
    pktab = with_pk ? pk_table(pk; tmax = pk_tmax) : PKTable([0.0, 1.0], [0.0, 0.0])
    pbctab = with_pk ? pbc_table(Ng; tmax = pk_tmax) : PKTable([0.0, 1.0], [0.0, 0.0])
    return FilterSetup(cfg, F, WSGCache(F, nfull), nfull, N, rows, Ng, G, keep, K, kappa,
                       noise_pass_fraction(F), Dspk, H, boundary_end_bins(nfull, cfg.m), pk, pktab, pbctab)
end
