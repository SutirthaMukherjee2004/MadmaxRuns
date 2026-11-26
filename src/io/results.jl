function save_result(path::AbstractString, data)
    mkpath(dirname(abspath(path)))
    tmp = path * ".tmp"
    serialize(tmp, data); mv(tmp, path; force = true)
    return path
end
load_result(path::AbstractString) = deserialize(path)

function write_csv(path::AbstractString, header::AbstractVector, rows::AbstractMatrix)
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        println(io, join(header, ","))
        writedlm(io, rows, ',')
    end
    return path
end
