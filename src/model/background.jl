const F0_BKG = 4.218e6
const BKG_UNIT = 1e22

function background(r::AbstractVector{<:Real}; grid::SpectrumGrid = PAPER_GRID)
    length(r) == 8 || throw(ArgumentError("the background needs 8 parameters r_1..r_8"))
    f0 = F0_BKG
    b = Vector{Float64}(undef, grid.nbins)
    @inbounds for (i, f) in enumerate(frequencies(grid))
        beta = (f - 25e6 * (1 + r[1] / 15)) / (20e6 * (1 + r[2] / 10))
        v = 1e-20 * (erf((f - f0) / 5e6) * (f0 / f)^3 + exp(-beta^2)) +
            4e-22 * (1 + r[3]) * sin((f + r[4] * f0) / 2.5e6) +
            5e-24 * ((1 + r[5]) * sin((f + r[6] * f0) / 0.25e6) + (1 + r[7]) * sin((f + r[8] * f0) / 0.1e6))
        b[i] = v * BKG_UNIT
    end
    return b
end

background(rng::AbstractRNG; grid::SpectrumGrid = PAPER_GRID) = background(randn(rng, 8); grid)

nominal_background(; grid::SpectrumGrid = PAPER_GRID) = background(zeros(8); grid)

median_background(; grid::SpectrumGrid = PAPER_GRID) = median(nominal_background(; grid))
