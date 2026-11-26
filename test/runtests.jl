using Test, MadmaxRuns, LinearAlgebra, Random, Statistics
using MadmaxRuns: rfft, irfft
const FULL = get(ENV, "MADMAX_TEST_FULL", "0") == "1"

@testset "MadmaxRuns" begin
    include("test_model.jl")
    include("test_filter.jl")
    include("test_noise.jl")
    include("test_whitening.jl")
    include("test_boundary.jl")
    include("test_stats.jl")
    include("test_chain.jl")
    include("test_checkpoints.jl")
    include("test_studies.jl")
end
