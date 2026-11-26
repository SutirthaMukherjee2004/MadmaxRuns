using MadmaxRuns, MadmaxRunsPlots

paths = prepare_output_layout()
require_all = "--require-all" in ARGS
paper = "--paper" in ARGS
dirs = filter(arg -> !(arg in ("--require-all", "--paper")), ARGS)
length(dirs) <= 2 || throw(ArgumentError("usage: make_figures.jl [data_dir] [paperplots_dir] [--require-all] [--paper]"))
resultdir = isempty(dirs) ? paths.data : dirs[1]
figdir = length(dirs) == 2 ? dirs[2] : paths.paperplots
with_run_log(joinpath(paths.logs, "paperplots.log")) do
    log_message("Rendering studies from $resultdir to $figdir")
    made = paper ? make_paper(resultdir, figdir) : make_all(resultdir, figdir; require_all)
    isempty(made) && error("no paper figures were rendered")
    log_message("$(length(made)) figures written to $figdir")
end
