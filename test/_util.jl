using ITensors

function _random_mpo(
    rng::AbstractRNG, sites::AbstractVector{<:AbstractVector{Index{T}}}; m::Int=1
) where {T}
    sites_ = collect(Iterators.flatten(sites))
    Ψ = random_mps(rng, sites_, m)
    tensors = ITensor[]
    pos = 1
    for i in 1:length(sites)
        push!(tensors, prod(Ψ[pos:(pos + length(sites[i]) - 1)]))
        pos += length(sites[i])
    end
    return MPO(tensors)
end

function _random_mpo(sites::AbstractVector{<:AbstractVector{Index{T}}}; m::Int=1) where {T}
    return _random_mpo(Random.default_rng(), sites; m=m)
end
