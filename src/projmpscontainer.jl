
# FIXME: do you need this struct? We can just use a vector of ProjMPS objects.
struct ProjMPSContainer
    data::Vector{ProjMPS}
    sites::Vector{Vector{Index}}
    projector::Projector

    function ProjMPSContainer(data::AbstractVector{ProjMPS})
        for n in 2:length(data)
            data[n].sites == data[1].sites ||
                error("Sitedims mismatch $(data[n].sites) != $(data[1].sites)")
        end
        projector = reduce(|, x.projector for x in data)
        return new(data, data[1].sites, projector)
    end
end

function Quantics.extractdiagonal(projmpss::ProjMPSContainer, sites)
    return ProjMPSContainer([
        Quantics.extractdiagonal(projmps, sites) for projmps in projmpss.data
    ])
end

function Quantics.rearrange_siteinds(projmpss::ProjMPSContainer, sites)
    return ProjMPSContainer([
        Quantics.rearrange_siteinds(projmps, sites) for projmps in projmpss.data
    ])
end

function ITensors.prime(Ψ::ProjMPSContainer, args...; kwargs...)
    return ProjMPSContainer([prime(projmps, args...; kwargs...) for projmps in Ψ.data])
end

function Quantics.makesitediagonal(projmpss::ProjMPSContainer, sites)
    return ProjMPSContainer([
        Quantics.makesitediagonal(projmps, sites) for projmps in projmpss.data
    ])
end
