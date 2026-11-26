using MadmaxRuns, Test
const ROOT = dirname(@__DIR__)
outdir = isempty(ARGS) ? joinpath(ROOT, "smoke") : ARGS[1]

include(joinpath(ROOT, "test", "runtests.jl"))

failed = String[]
for name in STUDY_ORDER
    t = time()
    try
        run_study(name; quick = true, outdir)
        println("smoke: $name ok (", round(time() - t, digits = 1), " s)")
    catch e
        push!(failed, name)
        println("smoke: $name FAILED"); showerror(stdout, e, catch_backtrace()); println()
    end
end
if isempty(failed)
    println("ALL TESTS PASSED")
else
    println("FAILED STUDIES: ", join(failed, ", ")); exit(1)
end
