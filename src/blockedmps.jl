
struct BlockedMPS
    data::OrderedDict{Projector,ProjMPS}

    function BlockedMPS(data::AbstractVector{ProjMPS})
        sites_all = [siteinds(prjmps) for prjmps in data]
        for n in 2:length(data)
            Set(sites_all[n]) == Set(sites_all[1]) || error("Sitedims mismatch")
        end
        for (n, a) in enumerate(data), (m, b) in enumerate(data)
            if n != m
                if hasoverlap(a.projector, b.projector)
                    error("$(a.projector) and $(b.projector) is overlapping")
                end
            end
        end
        dict_ = OrderedDict{Projector,ProjMPS}(
            data[i].projector => data[i] for i in 1:length(data)
        )
        return new(dict_)
    end
end

BlockedMPS(data::ProjMPS) = BlockedMPS([data])

function ITensors.siteinds(obj::BlockedMPS)
    return [Set(x) for x in ITensors.siteinds(first(values(obj.data)))]
end

Base.length(obj::BlockedMPS) = length(obj.data)

Base.getindex(bmps::BlockedMPS, i::Integer)::ProjMPS =
    first(Iterators.drop(values(bmps.data), i - 1))

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
    return BlockedMPS([extractdiagonal(prjmps, site) for prjmps in values(obj)])
end

function Quantics.rearrange_siteinds(obj::BlockedMPS, sites)
    return BlockedMPS([
        Quantics.rearrange_siteinds(prjmps, sites) for prjmps in values(obj)
    ])
end

function ITensors.prime(Ψ::BlockedMPS, args...; kwargs...)
    return BlockedMPS([prime(prjmps, args...; kwargs...) for prjmps in values(Ψ.data)])
end

function ITensors.norm(M::BlockedMPS)
    return sqrt(reduce(+, (x^2 for x in ITensors.norm.(values(M)))))
end

function Quantics.makesitediagonal(obj::BlockedMPS, site)
    return BlockedMPS([
        _makesitediagonal(prjmps, site; baseplev=baseplev) for prjmps in values(obj)
    ])
end

function _makesitediagonal(obj::BlockedMPS, site; baseplev=0)
    return BlockedMPS([
        _makesitediagonal(prjmps, site; baseplev=baseplev) for prjmps in values(obj)
    ])
end

Base.getindex(obj::BlockedMPS, p::Projector) = obj.data[p]

"""
Add two BlockedMPS objects.

If the two projects have the same projectors in the same order, the resulting BlockedMPS will have the same projectors in the same order.
"""
function Base.:+(
    a::BlockedMPS,
    b::BlockedMPS;
    alg=ITensors.Algorithm"directsum"(),
    cutoff=0.0,
    maxdim=typemax(Int),
    kwargs...,
)::BlockedMPS
    alg = ITensors.Algorithm(alg)
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

function ITensors.truncate(obj::BlockedMPS; kwargs...)::BlockedMPS
    return BlockedMPS([truncate(v; kwargs...) for v in values(obj)])
end

# Only for debug
function ITensors.MPS(obj::BlockedMPS; cutoff=1e-25, maxdim=typemax(Int), kwargs...)::MPS
    return reduce((x, y) -> +(x, y; kwargs), values(obj.data)).data # direct sum
end
