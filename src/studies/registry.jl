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

function study_options(kwargs)
    fields = fieldnames(StudyOptions)
    okw = Dict{Symbol,Any}(k => v for (k, v) in kwargs if k in fields)
    skw = Dict{Symbol,Any}(k => v for (k, v) in kwargs if !(k in fields))
    o = StudyOptions(; okw...)
    o.quick && quick!(o)
    return o, skw
end

function study_slot(name::AbstractString, o::StudyOptions)
    name == "all" && return "all"
    haskey(STUDIES, name) || throw(ArgumentError("unknown study '$name'; available: " * join(sort(collect(keys(STUDIES))), ", ")))
    if name in ("table_a", "table_b")
        foreach(parse_shape, o.shapes)
        return name * shape_suffix(o)
    end
    return name
end

function study_key(name, o, skw)
    io = IOBuffer()
    options = Tuple(getfield(o, field) for field in fieldnames(StudyOptions) if field != :outdir)
    extra_options = sort!(collect(pairs(skw)); by=first)
    serialize(io, (1, workflow_fingerprint(), name, options, extra_options))
    return bytes2hex(SHA.sha256(take!(io)))
end

study_file_digest(path) = open(io -> bytes2hex(SHA.sha256(io)), path)

function read_study_receipt(receipt, outputs, key)
    isfile(receipt) && all(isfile, outputs) || return nothing
    try
        data = load_result(receipt)
        data.key == key && data.digests == study_file_digest.(outputs) || return nothing
        return load_result(first(outputs))
    catch err
        err isa InterruptException && rethrow()
        log_message("checkpoint $receipt could not be read: $(sprint(showerror, err)); recomputing study")
        return nothing
    end
end

function run_study(name::AbstractString; resume::Bool = false, kwargs...)
    o, skw = study_options(kwargs)
    slot = study_slot(name, o)
    mkpath(o.outdir)
    if name == "all"
        return Dict(n => run_study(n; resume, kwargs...) for n in STUDY_ORDER)
    end
    if o.quick
        for (k, v) in get(QUICK_KW, name, Dict{Symbol,Any}())
            haskey(skw, k) || (skw[k] = v)
        end
    end
    checkpoint_dir = joinpath(o.outdir, ".study")
    mkpath(checkpoint_dir)
    receipt = joinpath(checkpoint_dir, slot * ".jls")
    outputs = [result_path(o, slot)]
    name in ("table_a", "table_b") && push!(outputs, joinpath(o.outdir, slot * ".csv"))
    try
        return with_checkpoint_lock(joinpath(checkpoint_dir, slot * ".lock")) do
            if name in ("correlation", "pvalues")
                o.tau_min, o.tau0_homo, o.tau0_hetero = thresholds(o)
            elseif name in ("fptp", "tau0")
                o.tau_min = minimum_threshold(o)
            end
            key = study_key(name, o, skw)
            if resume
                result = read_study_receipt(receipt, outputs, key)
                if result !== nothing
                    log_message("study $slot resumed from verified result $(first(outputs))")
                    return result
                end
            end
            log_message("study $slot: n = $(o.n), m = $(o.m), $(o.kind), nspec = $(o.nspec), nnull = $(o.nnull), nfilt = $(o.nfilt), threads = $(Threads.nthreads())")
            t = time()
            result = if name == "pvalues"
                progress = CheckpointRun(joinpath(checkpoint_dir, "pvalues"), (; key); resume)
                log_message("p-value study checkpoints: $(progress.dir)")
                study_pvalues(o; skw..., progress)
            else
                STUDIES[name](o; skw...)
            end
            save_result(receipt, (; key, digests = study_file_digest.(outputs)))
            log_message(@sprintf("study %s finished in %.1f s \u2192 %s", slot, time() - t, o.outdir))
            return result
        end
    catch err
        log_message("study $slot stopped: $(sprint(showerror, err, catch_backtrace()))")
        rethrow()
    end
end
