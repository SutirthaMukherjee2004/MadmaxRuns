Base.@kwdef struct PipelineConfig
    n::Int = 4
    m::Int = 64
    kind::Symbol = :SG
    tau::Float64 = 1e-7
    drop_nyquist::Bool = true
    W::Int = 96
    estimator::Symbol = :window
    tol::Float64 = 1e-10
    maxit::Int = 60
    boundary_removal::Bool = true
    nboundary::Int = 6
    pk_order::Int = 2
    qmax::Int = 60
    alpha::Float64 = 0.05
    nblock::Int = 100
    grid::SpectrumGrid = PAPER_GRID
end
with_tau(c::PipelineConfig, tau::Real) = PipelineConfig(; (k => getfield(c, k) for k in fieldnames(PipelineConfig))..., tau = Float64(tau))
function with(c::PipelineConfig; kw...)
    d = Dict{Symbol,Any}(k => getfield(c, k) for k in fieldnames(PipelineConfig))
    for (k, v) in kw
        d[k] = v
    end
    return PipelineConfig(; d...)
end
