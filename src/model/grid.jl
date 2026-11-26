struct SpectrumGrid
    nbins::Int
    df_hz::Float64
    fref_hz::Float64
end

const PAPER_GRID = SpectrumGrid(25_000, 2e3, 4.218e6)

frequencies(g::SpectrumGrid) = g.fref_hz .+ g.df_hz .* (0:g.nbins-1)

frequency_bin(g::SpectrumGrid, offset_hz::Real) = round(Int, offset_hz / g.df_hz) + 1
