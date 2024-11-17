
"""
BlockedMPS is a structure that holds multiple MPSs (ProjMPS) that are associated with different non-overlapping projectors.
"""
struct BlockedMPS
    data::OrderedDict{Projector,ProjMPS}

    function BlockedMPS(data::AbstractVector{ProjMPS})
        sites_all = [siteinds(prjmps) for prjmps in data]
        for n in 2:length(data)
            Set(sites_all[n]) == Set(sites_all[1]) || error("Sitedims mismatch")
        end
        isdisjoint([prjmps.projector for prjmps in data]) || error("Projectors are overlapping")

        dict_ = OrderedDict{Projector,ProjMPS}(
            data[i].projector => data[i] for i in 1:length(data)
        )
        return new(dict_)
    end
end

BlockedMPS(data::ProjMPS) = BlockedMPS([data])

"""
Return the site indices of the BlockedMPS.
The site indices are returned as a vector of sets, where each set corresponds to the site indices at each site.
"""
function siteindices(obj::BlockedMPS)
    return [Set(x) for x in ITensors.siteinds(first(values(obj.data)))]
end

ITensors.siteinds(obj::BlockedMPS) = siteindices(obj)

"""
Get the number of sites in the BlockedMPS
"""
Base.length(obj::BlockedMPS) = length(first(obj.data))

"""
Indexing for BlockedMPS. This is deprecated and will be removed in the future.
"""
function Base.getindex(bmps::BlockedMPS, i::Integer)::ProjMPS
    @warn "Indexing for BlockedMPS is deprecated. Use getindex(bmps, p::Projector) instead."
    return first(Iterators.drop(values(bmps.data), i - 1))
end

Base.getindex(obj::BlockedMPS, p::Projector) = obj.data[p]

function Base.iterate(bmps::BlockedMPS, state)
    return iterate(bmps.data, state)
end

function Base.iterate(bmps::BlockedMPS)
    return iterate(bmps.data)
end

"""
Return the keys, i.e., projectors of the BlockedMPS.
"""
function Base.keys(obj::BlockedMPS)
    return keys(obj.data)
end

"""
Return the values, i.e., ProjMPS of the BlockedMPS.
"""
function Base.values(obj::BlockedMPS)
    return values(obj.data)
end

"""
Extract diagonal of the BlockedMPS for `s`, `s'`, ... for a given site index `s`,
where `s` must have a prime level of 0.
"""
function extractdiagonal(obj::BlockedMPS, site)
    return BlockedMPS([extractdiagonal(prjmps, site) for prjmps in values(obj)])
end

"""
Rearrange the site indices of the BlockedMPS according to the given order.
If nessecary, tensors are fused or split to match the new order.
"""
function rearrange_siteinds(obj::BlockedMPS, sites)
    return BlockedMPS([rearrange_siteinds(prjmps, sites) for prjmps in values(obj)])
end

function ITensors.prime(Ψ::BlockedMPS, args...; kwargs...)
    return BlockedMPS([prime(prjmps, args...; kwargs...) for prjmps in values(Ψ.data)])
end

"""
Return the norm of the BlockedMPS.
"""
function ITensors.norm(M::BlockedMPS)
    return sqrt(reduce(+, (x^2 for x in ITensors.norm.(values(M)))))
end

"""
Make the BlockedMPS diagonal for a given site index `s` by introducing a dummy index `s'`.
"""
function makesitediagonal(obj::BlockedMPS, site)
    return BlockedMPS([
        _makesitediagonal(prjmps, site; baseplev=baseplev) for prjmps in values(obj)
    ])
end

function _makesitediagonal(obj::BlockedMPS, site; baseplev=0)
    return BlockedMPS([
        _makesitediagonal(prjmps, site; baseplev=baseplev) for prjmps in values(obj)
    ])
end

"""
Add two BlockedMPS objects.

If the two projects have the same projectors in the same order, the resulting BlockedMPS will have the same projectors in the same order.
By default, we use `directsum` algorithm to compute the sum and no truncation is performed.
"""
function Base.:+(
    a::BlockedMPS,
    b::BlockedMPS;
    alg="directsum",
    cutoff=0.0,
    maxdim=typemax(Int),
    kwargs...,
)::BlockedMPS
    data = ProjMPS[]
    for k in unique(vcat(collect(keys(a)), collect(keys(b)))) # preserve order
        if k ∈ keys(a) && k ∈ keys(b)
            a[k].projector == b[k].projector || error("Projectors mismatch at $(k)")
            push!(data, +(a[k], b[k]; alg, cutoff, maxdim, kwargs...))
        elseif k ∈ keys(a)
            push!(data, a[k])
        elseif k ∈ keys(b)
            push!(data, b[k])
        else
            error("Something went wrong")
        end
    end
    return BlockedMPS(data)
end

function Base.:*(a::BlockedMPS, b::Number)::BlockedMPS
    return BlockedMPS([a[k] * b for k in keys(a)])
end

function Base.:*(a::Number, b::BlockedMPS)::BlockedMPS
    return b * a
end

function Base.:-(obj::BlockedMPS)::BlockedMPS
    return -1 * obj
end

function truncate(obj::BlockedMPS; kwargs...)::BlockedMPS
    return BlockedMPS([truncate(v; kwargs...) for v in values(obj)])
end

# Only for debug
function ITensorMPS.MPS(obj::BlockedMPS; cutoff=1e-25, maxdim=typemax(Int))::MPS
    return reduce(
        (x, y) -> truncate(+(x, y; alg="directsum"); cutoff, maxdim), values(obj.data)
    ).data # direct sum
end

function ITensorMPS.MPO(obj::BlockedMPS; cutoff=1e-25, maxdim=typemax(Int))::MPO
    return MPO(collect(MPS(obj; cutoff=cutoff, maxdim=maxdim, kwargs...)))
end
