using MadmaxRuns, JSON, Base64

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const PAPER_NOTEBOOK = joinpath(ROOT, "notebooks", "paperplot.ipynb")

function run_paper_notebook(path::AbstractString = PAPER_NOTEBOOK; template::AbstractString = PAPER_NOTEBOOK)
    notebook = JSON.parsefile(template)
    cells = notebook["cells"]
    length(cells) == 22 && all(c -> c["cell_type"] == "code", cells) ||
        throw(ArgumentError("The paper notebook must have one setup cell and 21 figure cells"))
    execution = Module(:PaperNotebook)
    cd(ROOT) do
        for (i, cell) in enumerate(cells)
            value = Base.include_string(execution, join(cell["source"]), "$(basename(template)):cell_$i")
            cell["execution_count"] = i
            cell["outputs"] = if i == 1
                Any[]
            else
                mime = MIME"image/png"()
                Base.invokelatest(showable, mime, value) ||
                    error("Figure $(i - 1) did not return a PNG; the notebook has not been updated")
                png = Base.invokelatest(repr, mime, value)
                log_message("Rendered and embedded paper Figure $(i - 1)")
                [Dict("output_type" => "display_data", "metadata" => Dict(),
                      "data" => Dict("image/png" => base64encode(png)))]
            end
        end
    end
    MadmaxRuns.atomic_write(path) do io
        JSON.print(io, notebook, 1)
        println(io)
    end
    return path
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) <= 1 || throw(ArgumentError("usage: run_paper_notebook.jl [output_notebook]"))
    output = isempty(ARGS) ? PAPER_NOTEBOOK : abspath(only(ARGS))
    paths = prepare_output_layout()
    with_run_log(joinpath(paths.logs, "paperplots.log")) do
        log_message("Building paper notebook from $(paths.data)")
        run_paper_notebook(output)
        log_message("Saved $(output)")
    end
end
