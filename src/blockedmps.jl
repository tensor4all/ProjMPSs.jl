
struct BlockedMPS{K}
    data::Dict{K,ProjMPS}

    function BlockedMPS{K}(data::Dict{K,ProjMPS})::BlockedMPS{K} where {K}
        sites_all = [siteinds(prjmps) for (_, prjmps) in data]
        for n in 2:length(data)
            Set(sites_all[n]) == Set(sites_all[1]) || error("Sitedims mismatch")
        end
        projectors = [prjmps.projector for (_, prjmps) in data]
        for (n, a) in enumerate(projectors), (m, b) in enumerate(projectors)
            if n != m
                if hasoverlap(a, b)
                    error("$(a) and $(b) is overlapping")
                end
            end
        end
        return new{K}(data)
    end

end

function BlockedMPS(data::Dict{K,ProjMPS})::BlockedMPS{K} where {K}
    return BlockedMPS{K}(data)
end

function BlockedMPS(data::AbstractVector{ProjMPS})::BlockedMPS
    dict = Dict(data[i].projector => data[i] for i in 1:length(data))
    if length(dict) != length(data)
        error("Projectors are not unique")
    end
    return BlockedMPS(dict)
end

BlockedMPS(data::ProjMPS) = BlockedMPS([data])

ITensors.siteinds(obj::BlockedMPS) = ITensors.siteinds(first(obj.data))

Base.length(obj::BlockedMPS) = length(obj.data)

"""
Iterate over all data like a dictionary
"""
function Base.iterate(obj::BlockedMPS, state=1)
    if state > length(obj.data)
        return nothing
    end
    return (obj.data[state], state + 1)
end

function Base.keys(obj::BlockedMPS)
    return keys(obj.data)
end

function Base.values(obj::BlockedMPS)
    return values(obj.data)
end

function extractdiagonal(obj::BlockedMPS{K}, site) where {K}
    return BlockedMPS{K}(Dict(k => extractdiagonal(prjmps, site) for (k, prjmps) in obj))
end

function Quantics.rearrange_siteinds(obj::BlockedMPS{K}, sites) where {K}
    return BlockedMPS{K}(
        Dict(k => Quantics.rearrange_siteinds(projmps, sites) for (k, prjmps) in obj)
    )
end

function ITensors.prime(Ψ::BlockedMPS{K}, args...; kwargs...) where {K}
    return BlockedMPS{K}(Dict(k => prime(prjmps, args...; kwargs...) for (k, prjmps) in Ψ))
end

function Quantics.makesitediagonal(obj::BlockedMPS{K}, site) where {K}
    return BlockedMPS{K}(
        Dict(k => _makesitediagonal(prjmps, site; baseplev=baseplev) for (k, prjmps) in obj)
    )
end

function Quantics.makesitediagonal(obj::BlockedMPS{K}, site) where {K}
    return BlockedMPS{K}(Dict(k => Quantics.makesitediagonal(prjmps, site) for (k, prjmps) in obj))
end

Base.getindex(obj::BlockedMPS{K}, i::K) where {K} = obj.data[i]

#function Base.getindex(obj::BlockedMPS, p::Projector)
    #idx = findfirst(x -> x.projector == p, obj.data)
    #idx !== nothing || error("Projector $(p) not found")
    #return first(Iterators.drop(obj.data, idx - 1))
#end

function Base.:+(a::BlockedMPS{K}, b::BlockedMPS{K}; alg=ITensors.Algorithm"directsum"(), cutoff=0.0, maxdim=typemax(Int), kwargs...)::BlockedMPS{K} where {K}
    data = Dict{K,ProjMPS}()
    for k in union(keys(a), keys(b))
        if k ∈ keys(a) && k ∈ keys(b)
            a[k].projector == b[k].projector || error("Projectors mismatch at $(k)")
            data[k] = +(a[k], b[k]; cutoff=cutoff, maxdim=maxdim, kwargs...)
        elseif k ∈ keys(a)
            data[k] = a[k]
        elseif k ∈ keys(b)
            data[k] = b[k]
        else
            error("Something went wrong")
        end
    end
    return BlockedMPS{K}(data)
end

function Base.:*(a::BlockedMPS{K}, b::Number)::BlockedMPS{K} where {K}
    return BlockedMPS{K}(Dict(k => a[k] * b for k in keys(a)))
end

function Base.:*(a::Number, b::BlockedMPS)::BlockedMPS
    return b * a
end

function Base.:-(obj::BlockedMPS)::BlockedMPS
    return -1 * obj
end

function ITensors.truncate(obj::BlockedMPS; kwargs...)::BlockedMPS
    return BlockedMPS(Dict(k => truncate(v; kwargs...) for (k, v) in obj))
end

# Only for debug
function ITensors.MPS(obj::BlockedMPS; cutoff=1e-25, maxdim=typemax(Int), kwargs...)::MPS
    return reduce((x, y) -> +(x, y; kwargs), values(obj.data)).data # direct sum
end
