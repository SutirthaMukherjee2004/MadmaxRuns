using MadmaxRuns

function main(args)
    if isempty(args) || args[1] in ("-h", "--help", "list")
        println("studies: ", join(STUDY_ORDER, ", "), ", all")
        println("options: ", join(string.(fieldnames(StudyOptions)), ", "), ", resume, output_root, log")
        return
    end
    kw = parse_cli(args[2:end])
    paths = prepare_output_layout(pop!(kw, :output_root, get(ENV, "MADMAX_OUT", "results")))
    kw[:outdir] = get(kw, :outdir, paths.data)
    o, _ = MadmaxRuns.study_options(kw)
    slot = MadmaxRuns.study_slot(args[1], o)
    with_run_log(pop!(kw, :log, joinpath(paths.logs, "run_$slot.log"))) do
        log_message("Starting study $slot in $(paths.root)")
        configure_threads()
        run_study(args[1]; kw...)
    end
end

main(ARGS)
