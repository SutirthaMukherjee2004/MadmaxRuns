using MadmaxRunsPlots
resultdir = isempty(ARGS) ? "results" : ARGS[1]
figdir = length(ARGS) >= 2 ? ARGS[2] : joinpath(resultdir, "figs")
made = make_all(resultdir, figdir)
println(length(made), " figures written to ", figdir)
