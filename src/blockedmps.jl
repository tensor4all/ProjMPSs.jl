
struct BlockedMPS
    data::OrderedSet{ProjMPS}
    projectors::OrderedSet{Projector}

    function BlockedMPS(data::AbstractVector{ProjMPS})
        for n in 2:length(data)
            Set(siteinds(data[n])) == Set(siteinds(data[1])) || error("Sitedims mismatch")
        end
        for n in 1:length(data), m in 1:length(data)
            if n != m
                if hasoverlap(data[n].projector, data[m].projector)
                    error("$(data[n]) and $(data[m]) have overlapping projectors")
                end
            end
        end
        return new(OrderedSet(data), OrderedSet([proj.projector for proj in data]))
    end
end

BlockedMPS(data::ProjMPS) = BlockedMPS([data])

ITensors.siteinds(obj::BlockedMPS) = ITensors.siteinds(first(obj.data))

Base.length(obj::BlockedMPS) = length(obj.data)

function Base.iterate(obj::BlockedMPS, state=1)
    if state > length(obj.data)
        return nothing
    end
    return Base.iterate(obj.data, state)
    #return (obj.data[state], state + 1)
end

function Quantics.extractdiagonal(obj::BlockedMPS, site)
    return BlockedMPS([Quantics.extractdiagonal(projmps, site) for projmps in obj.data])
end

function Quantics.rearrange_siteinds(obj::BlockedMPS, sites)
    return BlockedMPS([Quantics.rearrange_siteinds(projmps, sites) for projmps in obj.data])
end

function ITensors.prime(Ψ::BlockedMPS, args...; kwargs...)
    return BlockedMPS([prime(projmps, args...; kwargs...) for projmps in Ψ.data])
end

function Quantics.makesitediagonal(obj::BlockedMPS, site; plev=0)
    return BlockedMPS([Quantics.makesitediagonal(projmps, site; plev=0) for projmps in obj.data])
end

Base.getindex(obj::BlockedMPS, i::Integer) = obj.data[i]

function Base.getindex(obj::BlockedMPS, p::Projector)
    idx = findfirst(x -> x.projector == p, obj.data)
    idx !== nothing || error("Projector $(p) not found")
    return first(Iterators.drop(obj.data, idx - 1))
end

function Base.:+(a::BlockedMPS, b::BlockedMPS)::BlockedMPS
    projmpss = ProjMPS[]
    for p in union(a.projectors, b.projectors)
        if p ∈ a.projectors && p ∈ b.projectors
            push!(projmpss, a[p] + b[p])
        elseif p ∈ a.projectors
            push!(projmpss, a[p])
        elseif p ∈ b.projectors
            push!(projmpss, b[p])
        else
            error("Something went wrong")
        end
    end
    return BlockedMPS(projmpss)
end

function Base.:*(a::BlockedMPS, b::Number)::ProjMPS
    return BlockedMPS([a[p] * b for p in a.projector])
end

function Base.:*(a::Number, b::BlockedMPS)::ProjMPS
    return BlockedMPS([a * b[p] for p in b.projector])
end

function Base.:-(obj::BlockedMPS)::ProjMPS
    return -1 * obj
end

function ITensors.truncate(obj::BlockedMPS; kwargs...)::ProjMPS
    return BlockedMPS([truncate(projmps; kwargs...) for projmps in obj.data])
end

# Only for debug
function ITensors.MPS(obj::BlockedMPS)
    return reduce(+, obj.data).data # direct sum
end
