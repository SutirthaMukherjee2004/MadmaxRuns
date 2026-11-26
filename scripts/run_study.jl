using MadmaxRuns

function main(args)
    if isempty(args) || args[1] in ("-h", "--help", "list")
        println("studies: ", join(STUDY_ORDER, ", "), ", all")
        println("options: ", join(string.(fieldnames(StudyOptions)), ", "))
        return
    end
    kw = parse_cli(args[2:end])
    outdir = get(kw, :outdir, "results"); mkpath(outdir)
    open_log(joinpath(outdir, "run_$(args[1]).log"))
    configure_threads()
    try
        run_study(args[1]; kw...)
    finally
        close_log()
    end
end

main(ARGS)
