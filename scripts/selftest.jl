using MadmaxRuns, Test
const ROOT = dirname(@__DIR__)
paths = prepare_output_layout(isempty(ARGS) ? joinpath(ROOT, "smoke") : ARGS[1])

with_run_log(joinpath(paths.logs, "smoke.log")) do
    log_message("Starting self-test in $(paths.root)")
    include(joinpath(ROOT, "test", "runtests.jl"))
    for name in ("filter", "sigmafit", "iteration")
        run_study(name; quick = true, resume = true, nspec = 2, outdir = paths.data)
    end
    log_message("ALL TESTS PASSED")
end
