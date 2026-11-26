function prepare_output_layout(root::AbstractString = get(ENV, "MADMAX_OUT", "results"))
    isempty(root) && throw(ArgumentError("output root must not be empty"))
    base = abspath(root)
    paths = (; root = base, logs = joinpath(base, "logs"), data = joinpath(base, "data"),
             paperplots = joinpath(base, "paperplots"))
    foreach(mkpath, (paths.logs, paths.data, paths.paperplots))
    return paths
end

function atomic_write(f, path::AbstractString)
    target = abspath(path)
    dir = dirname(target)
    mkpath(dir)
    tmp, io = mktemp(dir)
    try
        f(io)
        close(io)
        Base.Filesystem.rename(tmp, target)
    finally
        isopen(io) && close(io)
        isfile(tmp) && rm(tmp)
    end
    return path
end

save_result(path::AbstractString, data) = atomic_write(io -> serialize(io, data), path)
load_result(path::AbstractString) = deserialize(path)

function write_csv(path::AbstractString, header::AbstractVector, rows::AbstractMatrix)
    atomic_write(path) do io
        println(io, join(header, ","))
        writedlm(io, rows, ',')
    end
end
