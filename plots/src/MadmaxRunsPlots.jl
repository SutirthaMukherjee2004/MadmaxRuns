module MadmaxRunsPlots

using CairoMakie, Statistics, Printf
import CairoMakie: Makie
using MadmaxRuns

include("style.jl")
include("figures_model.jl")
include("figures_tests.jl")
include("figures_chain.jl")
include("paper.jl")

export make_all, paper_theme, PAPER_PLOT_STUDIES
export PaperRender, render_paper!, make_paper

const PAPER_PLOT_STUDIES = ("filter", "background", "spectra", "sigmafit", "iteration", "nobkg",
                            "tauscan", "tau0", "correlation", "fptp", "pvalues", "examples",
                            "leakage", "tolerance", "boundary", "tcut")

function make_all(resultdir::AbstractString, figdir::AbstractString; require_all::Bool = false)
    mkpath(figdir)
    missing = [name for name in PAPER_PLOT_STUDIES if !isfile(joinpath(resultdir, name * ".jls"))]
    if !isempty(missing)
        message = "Missing study results in $resultdir: " * join(missing, ", ")
        require_all && throw(ArgumentError(message))
        log_message(message)
    end
    require_all && serif_fonts() === nothing &&
        error("Paper export requires the DejaVu Serif font family")
    set_theme!(paper_theme())
    made = String[]
    load(name) = (p = joinpath(resultdir, name * ".jls"); isfile(p) ? load_result(p) : nothing)
    function make(f, name, args...)
        try
            push!(made, savefig(f(args...), figdir, name))
            log_message("Wrote $name (PDF and PNG)")
        catch err
            log_message("Figure $name failed: $(sprint(showerror, err))")
            rethrow()
        end
    end
    tau0 = load("tau0")
    t0 = tau0 === nothing ? nothing : (homo = tau0.cross["homo"].tau0, hetero = tau0.cross["hetero"].tau0)
    (r = load("filter")) !== nothing && make(fig_filter, "fig_filter", r)
    (r = load("iteration")) !== nothing && make(fig_iteration, "iter_convergence", r)
    (r = load("background")) !== nothing && make(fig_background, "fig_background", r)
    if (r = load("nobkg")) !== nothing
        make(fig_nobkg_data, "nobkg_fig1", r); make(fig_nobkg_pvalues, "nobkg_fig2", r)
        make(fig_shapes_pvalues, "fig_shapes_pvals", r)
    end
    (r = load("spectra")) !== nothing && make(fig_spectra, "res_spectra", r)
    (r = load("sigmafit")) !== nothing && make(fig_sigmafit, "deck_nobkg", r)
    if (r = load("correlation")) !== nothing
        make(fig_correlation, "corr_homo_hetero7_random_bg", r); make(fig_corrdiff, "corr_diff", r)
    end
    (r = load("tauscan")) !== nothing && make(fig_tauscan, "deck_tauscan", r, t0)
    (r = load("fptp")) !== nothing && make(fig_fptp, "deck_three_v2", r)
    (r = load("leakage")) !== nothing && make(fig_leakage, "res_leakmap", r)
    if (r = load("examples")) !== nothing
        make(fig_examples, @sprintf("deck_examples_tau%.0e", r.tau), r)
    end
    (r = load("tolerance")) !== nothing && make(fig_tolerance, "fig_tolerance", r)
    if tau0 !== nothing
        make(fig_tau0_map, "tau_safe_map", tau0); make(fig_tau0_crossing, "tau_safe_crossing", tau0)
    end
    if (r = load("boundary")) !== nothing
        make(fig_boundary, "explain_projection2", r); make(fig_boundary_svd, "boundary_singular_values", r)
    end
    if (r = load("tcut")) !== nothing
        make(fig_tcut, "tcut_block_v", r); make(fig_g2nd, "appg_second_order", r)
    end
    if (r = load("pvalues")) !== nothing
        for null in (:filt, :iid), tag in r.models
            make(fig_pvalues, "deck_pvalues_$(null)_$(tag)", r, null, tag)
        end
    end
    return made
end

end
