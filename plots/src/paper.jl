const PAPER_FIGURES = (
    ("filter", fig_filter),
    ("iteration", fig_iteration),
    ("background", fig_background),
    ("nobkg", fig_nobkg),
    ("nobkg", fig_shapes_pvalues),
    ("spectra", fig_spectra),
    ("sigmafit", fig_sigmafit),
    ("correlation", fig_correlation),
    ("correlation", fig_corrdiff),
    ("tauscan", fig_tauscan),
    ("fptp", fig_fptp),
    ("leakage", fig_leakage),
    ("", nothing),
    ("examples", fig_examples),
    ("tolerance", fig_tolerance),
    ("tau0", fig_tau0_map),
    ("tau0", fig_tau0_crossing),
    ("boundary", fig_boundary),
    ("tcut", fig_tcut),
    ("tcut", fig_g2nd),
    ("pvalues", fig_pvalues_paper),
)

struct PaperImage
    path::String
end

Base.show(io::IO, ::MIME"image/png", image::PaperImage) = write(io, read(image.path))

struct PaperRender
    resultdir::String
    figdir::String
    diagram::String
    theme::Theme
    results::Dict{String,Any}
end

function PaperRender(resultdir::AbstractString, figdir::AbstractString; require_full::Bool = true)
    required = unique([name for (name, _) in PAPER_FIGURES if !isempty(name)])
    missing = [name for name in required if !isfile(joinpath(resultdir, name * ".jls"))]
    isempty(missing) || throw(ArgumentError("Missing studies for the paper: " * join(missing, ", ")))
    diagram = normpath(joinpath(@__DIR__, "..", "..", "docs", "figure-13-analysis-pipeline.png"))
    isfile(diagram) || throw(ArgumentError("Missing Figure 13 diagram: $diagram"))
    serif_fonts() === nothing && error("Paper export requires the DejaVu Serif font family")
    mkpath(figdir)
    render = PaperRender(abspath(resultdir), abspath(figdir), diagram, paper_theme(), Dict{String,Any}())
    require_full && check_paper_counts(render)
    return render
end

study_result(render::PaperRender, name::String) = get!(render.results, name) do
    load_result(joinpath(render.resultdir, name * ".jls"))
end

function check_paper_counts(render::PaperRender)
    for name in ("examples", "fptp", "pvalues", "correlation", "tauscan")
        r = study_result(render, name)
        if hasproperty(r, :nspec)
            r.nspec == 5000 || throw(ArgumentError("$name requires 5000 spectra per case, found $(r.nspec)"))
        end
    end
    nobkg = study_result(render, "nobkg")
    length(nobkg.p_noise) == 5000 || throw(ArgumentError("Figures 4 and 5 require 5000 background trials"))
    all(d -> length(d.p) == 5000, values(nobkg.shapes)) ||
        throw(ArgumentError("Figures 4 and 5 require 5000 trials for each line shape"))
    study_result(render, "iteration").nspec == 100 ||
        throw(ArgumentError("Figure 2 requires 100 spectra per case"))
    study_result(render, "tcut").nmc == 500_000 ||
        throw(ArgumentError("Figures G.19 and G.20 require 500000 reference spectra"))
    return nothing
end

function render_paper!(render::PaperRender, number::Integer)
    study, builder = PAPER_FIGURES[number]
    name = @sprintf("figure-%02d", number)
    path = joinpath(render.figdir, name * ".png")
    if number == 13
        MadmaxRuns.atomic_write(path) do io
            write(io, read(render.diagram))
        end
    else
        data = study_result(render, study)
        with_theme(render.theme) do
            figure = if number == 10
                tau0 = study_result(render, "tau0")
                builder(data, (homo = tau0.cross["homo"].tau0, hetero = tau0.cross["hetero"].tau0))
            else
                builder(data)
            end
            savefig(figure, render.figdir, name)
        end
    end
    return PaperImage(path)
end

function make_paper(resultdir::AbstractString, figdir::AbstractString)
    render = PaperRender(resultdir, figdir)
    return map(eachindex(PAPER_FIGURES)) do number
        image = render_paper!(render, number)
        log_message("Wrote paper Figure $number")
        image.path
    end
end
