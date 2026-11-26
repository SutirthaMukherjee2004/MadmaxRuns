using Pkg
const ROOT = dirname(@__DIR__)

Pkg.activate(ROOT)
Pkg.add(["FFTW", "SpecialFunctions", "RunStatistics"])
Pkg.instantiate()
Pkg.precompile()

if "plots" in ARGS
    Pkg.activate(joinpath(ROOT, "plots"))
    Pkg.develop(path = ROOT)
    Pkg.add("CairoMakie")
    Pkg.instantiate()
    Pkg.precompile()
end
println("setup finished")
