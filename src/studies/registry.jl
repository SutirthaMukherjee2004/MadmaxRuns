const STUDIES = Dict{String,Function}(
    "filter" => study_filter,
    "nobkg" => study_nobkg,
    "table_a" => study_table_a,
    "table_b" => study_table_b,
    "background" => study_background,
    "spectra" => study_spectra,
    "sigmafit" => study_sigmafit,
    "iteration" => study_iteration,
    "tauscan" => study_tauscan,
    "tau0" => study_tau0,
    "correlation" => study_correlation,
    "fptp" => study_fptp,
    "pvalues" => study_pvalues,
    "examples" => study_examples,
    "leakage" => study_leakage,
    "tolerance" => study_tolerance,
    "boundary" => study_boundary,
    "tcut" => study_tcut,
)

const STUDY_ORDER = ["filter", "background", "spectra", "sigmafit", "iteration", "nobkg", "tauscan", "tau0",
                     "correlation", "fptp", "pvalues", "examples", "leakage", "tolerance", "boundary", "tcut",
                     "table_a", "table_b"]

function quick!(o::StudyOptions)
    o.nspec = min(o.nspec, 24); o.nnull = min(o.nnull, 60); o.nfilt = min(o.nfilt, 24)
    isempty(o.taus) && (o.taus = [1e-7, 1e-5, 1e-4, 1e-3, 0.01, 0.1])
    o.amplitudes == collect(0.0:4.0:28.0) && (o.amplitudes = [0.0, 12.0, 28.0])
    o.widths == collect(2.0:2.0:12.0) && (o.widths = [2.0, 8.0])
    return o
end

const QUICK_KW = Dict{String,Dict{Symbol,Any}}(
    "iteration" => Dict(:nsteps => 4),
    "leakage" => Dict(:ns => 2:2:4, :ms => [32, 64], :nbkg => 2),
    "tau0" => Dict(:ns => 2:2:4, :ms => [32, 64]),
    "tcut" => Dict(:ntau => [195], :nmc => 400, :Mexact => 400, :tgrid => collect(5.0:5.0:30.0), :qexact => 8),
    "tolerance" => Dict(:nnoise => 8),
)

function run_study(name::AbstractString; kwargs...)
    fields = fieldnames(StudyOptions)
    okw = Dict{Symbol,Any}(k => v for (k, v) in kwargs if k in fields)
    skw = Dict{Symbol,Any}(k => v for (k, v) in kwargs if !(k in fields))
    o = StudyOptions(; okw...)
    o.quick && quick!(o)
    mkpath(o.outdir)
    if name == "all"
        return Dict(n => run_study(n; kwargs...) for n in STUDY_ORDER)
    end
    haskey(STUDIES, name) || throw(ArgumentError("unknown study '$name'; available: " * join(sort(collect(keys(STUDIES))), ", ")))
    if o.quick
        for (k, v) in get(QUICK_KW, name, Dict{Symbol,Any}())
            haskey(skw, k) || (skw[k] = v)
        end
    end
    log_message("study $name: n = $(o.n), m = $(o.m), $(o.kind), nspec = $(o.nspec), nnull = $(o.nnull), nfilt = $(o.nfilt), threads = $(Threads.nthreads())")
    t = time()
    res = STUDIES[name](o; skw...)
    log_message(@sprintf("study %s finished in %.1f s \u2192 %s", name, time() - t, o.outdir))
    return res
end
