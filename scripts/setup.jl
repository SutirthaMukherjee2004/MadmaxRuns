using Pkg
const ROOT = dirname(@__DIR__)

all(arg -> arg in ("plots", "notebooks"), ARGS) ||
    throw(ArgumentError("usage: julia scripts/setup.jl [plots] [notebooks]"))

Pkg.activate(ROOT)
Pkg.instantiate()
Pkg.precompile()

if "plots" in ARGS || "notebooks" in ARGS
    cd(joinpath(ROOT, "plots")) do
        Pkg.activate(".")
        if !isfile("Manifest.toml") && !isfile("Manifest-v$(VERSION.major).$(VERSION.minor).toml")
            Pkg.develop(path = "..")
        end
        Pkg.instantiate()
        Pkg.precompile()
    end
end

if "notebooks" in ARGS
    Pkg.activate(joinpath(ROOT, "notebooks"))
    Pkg.instantiate()
    Pkg.precompile()
    @eval using IJulia
    Base.invokelatest(IJulia.installkernel, "MadmaxRuns", "--project=$(joinpath(ROOT, "notebooks"))";
                      specname = "madmax-runs", displayname = "MadmaxRuns (Julia)")
end
println("setup finished")
