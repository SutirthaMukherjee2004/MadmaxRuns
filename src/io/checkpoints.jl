const CHECKPOINT_FORMAT = 1

function with_checkpoint_lock(f, path)
    Sys.isunix() || return f()
    open(path, "a") do io
        ccall(:flock, Cint, (Cint, Cint), fd(io), 6) == 0 ||
            error("cannot lock $path; another process may be writing checkpoints for this run")
        return f()
    end
end

function workflow_fingerprint()
    root = dirname(dirname(realpath(@__DIR__)))
    paths = ["Project.toml"]
    isfile(joinpath(root, "Manifest.toml")) && push!(paths, "Manifest.toml")
    for dir in ("src", "scripts", "batch")
        for (parent, _, files) in walkdir(joinpath(root, dir))
            for file in files
                any(ext -> endswith(file, ext), (".jl", ".sh", ".txt")) &&
                    push!(paths, relpath(joinpath(parent, file), root))
            end
        end
    end
    sources = [(path, bytes2hex(sha256(read(joinpath(root, path))))) for path in sort!(paths)]
    return checkpoint_digest((VERSION, Sys.ARCH, Sys.KERNEL, sources,
                              pkgversion(FFTW), pkgversion(RunStatistics),
                              pkgversion(parentmodule(loggamma)), pkgversion(CodecZlib)))
end

function checkpoint_digest(data)
    io = IOBuffer()
    serialize(io, data)
    return bytes2hex(sha256(take!(io)))
end

mutable struct CheckpointRun
    dir::String
    signature::String
    reuse::Bool
    previous::UInt128
end

function CheckpointRun(root::AbstractString, inputs; resume::Bool = true)
    signature = checkpoint_digest((CHECKPOINT_FORMAT, workflow_fingerprint(), inputs))
    dir = joinpath(abspath(root), signature)
    mkpath(dir)
    return CheckpointRun(dir, signature, resume, UInt128(0))
end

function write_checkpoint(run::CheckpointRun, stage::Symbol, data; failed::Bool = false)
    path = joinpath(run.dir, string(stage, failed ? ".failed.jls.gz" : ".jls.gz"))
    id = rand(RandomDevice(), UInt128)
    atomic_write(path) do io
        stream = GzipCompressorStream(io)
        try
            serialize(stream, (; format = CHECKPOINT_FORMAT, signature = run.signature,
                               stage, previous = run.previous, id, data))
        finally
            close(stream)
        end
    end
    failed || (run.previous = id)
    return path
end

function load_checkpoint_file(path::AbstractString)
    bytes = open(path, "r") do io
        stream = GzipDecompressorStream(io)
        try
            read(stream)
        finally
            close(stream)
        end
    end
    return deserialize(IOBuffer(bytes))
end

function read_checkpoint(run::CheckpointRun, stage::Symbol, path::AbstractString)
    record = load_checkpoint_file(path)
    record.format == CHECKPOINT_FORMAT && record.signature == run.signature &&
        record.stage == stage && record.previous == run.previous ||
        error("checkpoint metadata does not match this run")
    run.previous = record.id
    return record.data
end

run_stage(f, ::Nothing, ::Symbol) = f()

function run_stage(f, run::CheckpointRun, stage::Symbol)
    path = joinpath(run.dir, string(stage, ".jls.gz"))
    if run.reuse && isfile(path)
        try
            data = read_checkpoint(run, stage, path)
            log_message("Resumed stage $stage")
            return data
        catch err
            err isa InterruptException && rethrow()
            log_message("Cannot read checkpoint for $stage: $(sprint(showerror, err)); recomputing from here")
        end
    end
    run.reuse = false
    log_message("Starting stage $stage")
    try
        data = f()
        write_checkpoint(run, stage, data)
        log_message("Checkpoint saved for $stage")
        return data
    catch err
        log_message("Stage $stage failed: $(sprint(showerror, err))")
        rethrow()
    end
end
