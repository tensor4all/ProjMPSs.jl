
struct BlockedMPS{IndsT}
    data::Dict{Projector{IndsT},ProjMPS}

    function BlockedMPS{IndsT}(data::Dict{Projector{IndsT},ProjMPS}) where {IndsT}
        sites_all = [siteinds(prjmps) for (_, prjmps) in data]
        for n in 2:length(data)
            Set(sites_all[n]) == Set(sites_all[1]) || error("Sitedims mismatch")
        end
        for (p, prjmps) in data
            if prjmps.projector != p
                error("Projector mismatch")
            end
        end
        projectors = [prjmps.projector for (_, prjmps) in data]
        for (n, a) in enumerate(projectors), (m, b) in enumerate(projectors)
            if n != m
                if hasoverlap(a, b)
                    error("$(a) and $(b) is overlapping")
                end
            end
        end
        return new(data)
    end
end

function BlockedMPS(data::Dict{Projector{IndsT},ProjMPS}) where {IndsT}
    return BlockedMPS{IndsT}(data)
end

function BlockedMPS(data::AbstractVector{ProjMPS})
    dict = Dict(data[i].projector => data[i] for i in 1:length(data))
    if length(dict) != length(data)
        error("Projectors are not unique")
    end
    return BlockedMPS(dict)
end

BlockedMPS(data::ProjMPS) = BlockedMPS([data])

ITensors.siteinds(obj::BlockedMPS) = ITensors.siteinds(first(values(obj.data)))

Base.length(obj::BlockedMPS) = length(obj.data)

function Base.iterate(bmps::BlockedMPS, state)
    return iterate(bmps.data, state)
end

function Base.iterate(bmps::BlockedMPS)
    return iterate(bmps.data)
end

function Base.keys(obj::BlockedMPS)
    return keys(obj.data)
end

function Base.values(obj::BlockedMPS)
    return values(obj.data)
end

function extractdiagonal(obj::BlockedMPS, site)
    return BlockedMPS(Dict(k => extractdiagonal(prjmps, site) for (k, prjmps) in obj))
end

function Quantics.rearrange_siteinds(obj::BlockedMPS, sites)
    return BlockedMPS(
        Dict(k => Quantics.rearrange_siteinds(prjmps, sites) for (k, prjmps) in obj)
    )
end

function ITensors.prime(Ψ::BlockedMPS, args...; kwargs...)
    return BlockedMPS(Dict(k => prime(prjmps, args...; kwargs...) for (k, prjmps) in Ψ))
end

function Quantics.makesitediagonal(obj::BlockedMPS, site)
    return BlockedMPS(
        Dict(k => _makesitediagonal(prjmps, site; baseplev=baseplev) for (k, prjmps) in obj)
    )
end

function _makesitediagonal(obj::BlockedMPS, site; baseplev=0)
    return BlockedMPS(
        Dict(k => _makesitediagonal(prjmps, site; baseplev=baseplev) for (k, prjmps) in obj)
    )
end

function Base.getindex(obj::BlockedMPS{IndsT}, p::Projector{IndsT}) where {IndsT}
    return obj.data[p]
end

function Base.:+(
    a::BlockedMPS{IndsT},
    b::BlockedMPS{IndsT};
    alg=ITensors.Algorithm"directsum"(),
    cutoff=0.0,
    maxdim=typemax(Int),
    kwargs...,
)::BlockedMPS{IndsT} where {IndsT}
    alg = ITensors.Algorithm(alg)
    data = Dict{Projector{IndsT},ProjMPS}()
    for k in union(keys(a), keys(b))
        if k ∈ keys(a) && k ∈ keys(b)
            a[k].projector == b[k].projector || error("Projectors mismatch at $(k)")
            data[k] = +(a[k], b[k]; alg, cutoff, maxdim, kwargs...)
        elseif k ∈ keys(a)
            data[k] = a[k]
        elseif k ∈ keys(b)
            data[k] = b[k]
        else
            error("Something went wrong")
        end
    end
    return BlockedMPS{IndsT}(data)
end

function Base.:*(a::BlockedMPS{IndsT}, b::Number)::BlockedMPS{IndsT} where {IndsT}
    return BlockedMPS([a[k] * b for k in keys(a)])
end

function Base.:*(a::Number, b::BlockedMPS{IndsT})::BlockedMPS{IndsT} where {IndsT}
    return b * a
end

function Base.:-(obj::BlockedMPS{IndsT})::BlockedMPS{IndsT} where {IndsT}
    return -1 * obj
end

function ITensors.truncate(
    obj::BlockedMPS{IndsT}; kwargs...
)::BlockedMPS{IndsT} where {IndsT}
    return BlockedMPS{IndsT}(Dict(k => truncate(v; kwargs...) for (k, v) in obj))
end

# Only for debug
function ITensors.MPS(obj::BlockedMPS; cutoff=1e-25, maxdim=typemax(Int), kwargs...)::MPS
    return reduce((x, y) -> +(x, y; kwargs), values(obj.data)).data # direct sum
end
