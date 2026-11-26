@enum NoiseModel Homoscedastic Heteroscedastic

function parse_noise_model(s::AbstractString)
    t = lowercase(s)
    startswith(t, "homo") && return Homoscedastic
    startswith(t, "hetero") && return Heteroscedastic
    throw(ArgumentError("unknown noise model \"$s\" (use homo or hetero)"))
end
model_tag(m::NoiseModel) = m == Homoscedastic ? "homo" : "hetero"

function noise_parameters(model::NoiseModel, level::Real, bmed::Real)
    sigmamed = level * bmed
    model == Homoscedastic && return (sigmamed, 0.0)
    sigma0 = sigmamed / sqrt(5)
    return (sigma0, (sigmamed - sigma0) / bmed)
end

noise_sigma(sigma0::Real, sigma1::Real, b::AbstractVector) = sigma0 .+ sigma1 .* b
