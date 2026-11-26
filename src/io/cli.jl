function parse_value(T::Type, s::AbstractString)
    T <: AbstractVector{NoiseModel} && return [parse_noise_model(strip(v)) for v in split(s, ',')]
    T <: AbstractVector{<:AbstractString} && return [String(strip(v)) for v in split(s, ',')]
    T <: AbstractVector && return [parse(eltype(T), strip(v)) for v in split(s, ',')]
    if T === Bool
        lowercase(s) in ("1", "true", "yes") && return true
        lowercase(s) in ("0", "false", "no") && return false
        throw(ArgumentError("invalid boolean '$s'"))
    end
    T === Symbol && return Symbol(s)
    T <: AbstractString && return String(s)
    return parse(T, s)
end

function guess_value(s::AbstractString)
    occursin(',', s) && return [guess_value(strip(v)) for v in split(s, ',')]
    occursin(':', s) && (p = parse.(Int, split(s, ':')); return length(p) == 2 ? (p[1]:p[2]) : (p[1]:p[2]:p[3]))
    v = tryparse(Int, s); v !== nothing && return v
    v = tryparse(Float64, s); v !== nothing && return v
    lowercase(s) in ("true", "false") && return lowercase(s) == "true"
    return String(s)
end

function parse_cli(args::AbstractVector{<:AbstractString})
    kw = Dict{Symbol,Any}()
    types = Dict(zip(fieldnames(StudyOptions), fieldtypes(StudyOptions)))
    types[:resume] = Bool
    types[:output_root] = String
    types[:log] = String
    for a in args
        occursin('=', a) || throw(ArgumentError("argument '$a' is not of the form key=value"))
        k, v = split(a, '='; limit = 2); key = Symbol(strip(k))
        kw[key] = haskey(types, key) ? parse_value(types[key], strip(v)) : guess_value(strip(v))
    end
    return kw
end
