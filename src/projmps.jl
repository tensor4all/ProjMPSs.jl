"""
An MPS with a projector.
"""
struct ProjMPS
    data::MPS
    projector::Projector

    function ProjMPS(data::AbstractMPS, projector::Projector{IndsT}) where {IndsT}
        _iscompatible(projector, data) || error(
            "Incompatible projector and data. Even small numerical noise can cause this error.",
        )
        projector = _trim_projector(data, projector)
        return new(MPS([x for x in data]), projector)
    end
end

ITensors.siteinds(obj::ProjMPS) = siteinds(MPO([x for x in obj.data]))

_allsites(Ψ::AbstractMPS) = collect(Iterators.flatten(siteinds(MPO(collect(Ψ)))))
_allsites(Ψ::ProjMPS) = _allsites(Ψ.data)

function _trim_projector(obj::AbstractMPS, projector)
    sites = Set(_allsites(obj))
    newprj = deepcopy(projector)
    for (k, v) in newprj.data
        if !(k in sites)
            delete!(newprj.data, k)
        end
    end
    return newprj
end

function ProjMPS(Ψ::AbstractMPS)
    IndsT = typeof(first(siteinds(MPS(collect(Ψ)))))
    return ProjMPS(Ψ, Projector{IndsT}())
end

# Conversion Functions
ITensors.MPS(projΨ::ProjMPS) = projΨ.data

function project(tensor::ITensor, projector::Projector{IndsT}) where {IndsT}
    slice = Union{Int,Colon}[
        isprojectedat(projector, idx) ? projector[idx] : Colon() for idx in inds(tensor)
    ]
    data_org = Array(tensor, inds(tensor)...)
    data_trim = zero(data_org)
    data_trim[slice...] .= data_org[slice...]
    return ITensor(data_trim, inds(tensor)...)
end

function project(
    projΨ::ProjMPS, projector::Projector{IndsT}
)::Union{Nothing,ProjMPS} where {IndsT}
    newprj = projector & projΨ.projector
    if newprj === nothing
        return nothing
    end

    return ProjMPS(
        MPS([project(projΨ.data[n], newprj) for n in 1:length(projΨ.data)]), newprj
    )
end

function project(
    Ψ::AbstractMPS, projector::Projector{IndsT}
)::Union{Nothing,ProjMPS} where {IndsT}
    return project(ProjMPS(Ψ), projector)
end

function project(
    projΨ::ProjMPS, projector::Dict{IndsT,Int}
)::Union{Nothing,ProjMPS} where {IndsT}
    return project(projΨ, Projector(projector))
end

function _iscompatible(projector::Projector, tensor::ITensor)
    # Lazy impmentation
    return norm(project(tensor, projector) - tensor) == 0.0
end

function _iscompatible(projector::Projector, Ψ::AbstractMPS)
    return all((_iscompatible(projector, x) for x in Ψ))
end

function _makesitediagonal(
    projmps::ProjMPS, sites::AbstractVector{Index{IndsT}}; baseplev=0
) where {IndsT}
    M_ = deepcopy(MPO(collect(MPS(projmps))))
    for site in sites
        target_site::Int = only(findsites(M_, site))
        M_[target_site] = _asdiagonal(M_[target_site], site; baseplev=baseplev)
    end
    return project(M_, projmps.projector)
end

function _makesitediagonal(projmps::ProjMPS, site::Index; baseplev=0)
    return _makesitediagonal(projmps, [site]; baseplev=baseplev)
end

# Quantics Functions
function Quantics.makesitediagonal(projmps::ProjMPS, site::Index)
    return _makesitediagonal(projmps, site; baseplev=0)
end

function Quantics.makesitediagonal(projmps::ProjMPS, sites::AbstractVector{Index})
    return _makesitediagonal(projmps, sites; baseplev=0)
end

function Quantics.makesitediagonal(projmps::ProjMPS, tag::String)
    mps_diagonal = Quantics.makesitediagonal(MPS(projmps), tag)
    projmps_diagonal = ProjMPS(mps_diagonal)

    target_sites = Quantics.findallsiteinds_by_tag(
        unique(noprime.(Iterators.flatten(siteinds(projmps)))); tag=tag
    )

    newproj = deepcopy(projmps.projector)
    for s in target_sites
        if isprojectedat(projmps.projector, s)
            newproj[s'] = newproj[s]
        end
    end

    return project(projmps_diagonal, newproj)
end

# FIXME: may be type unstable
function _find_site_allplevs(tensor::ITensor, site::Index; maxplev=10)
    plev(site) == 0 || error("Site index must be unprimed.")
    return [prime(site, plev) for plev in 0:maxplev if prime(site, plev) ∈ inds(tensor)]
end

function extractdiagonal(
    projmps::ProjMPS, sites::AbstractVector{Index{IndsT}}
) where {IndsT}
    tensors = collect(projmps.data)
    for i in eachindex(tensors)
        for site in intersect(sites, inds(tensors[i]))
            sitewithallplevs = _find_site_allplevs(tensors[i], site)
            tensors[i] = if length(sitewithallplevs) > 1
                tensors[i] = Quantics._extract_diagonal(tensors[i], sitewithallplevs...)
            else
                tensors[i]
            end
        end
    end

    projector = deepcopy(projmps.projector)
    for site in sites
        if site' in keys(projector.data)
            delete!(projector.data, site')
        end
    end
    return ProjMPS(MPS(tensors), projector)
end

function extractdiagonal(projmps::ProjMPS, site::Index{IndsT}) where {IndsT}
    return Quantics.extractdiagonal(projmps, [site])
end

function Quantics.extractdiagonal(projmps::ProjMPS, tag::String)::ProjMPS
    targetsites = Quantics.findallsiteinds_by_tag(
        unique(noprime.(_allsites(projmps))); tag=tag
    )
    return extractdiagonal(projmps, targetsites)
end

function Quantics.rearrange_siteinds(projmps::ProjMPS, sites)
    mps_rearranged = Quantics.rearrange_siteinds(MPS(projmps), sites)
    return project(ProjMPS(mps_rearranged), projmps.projector)
end

# Miscellaneous Functions
function Base.show(io::IO, obj::ProjMPS)
    return print(io, "ProjMPS projected on $(obj.projector.data)")
end

function ITensors.prime(Ψ::ProjMPS, args...; kwargs...)
    return ProjMPS(
        prime(MPS(Ψ), args...; kwargs...),
        prime.(siteinds(Ψ), args...; kwargs...),
        Ψ.projector,
    )
end

Base.isapprox(x::ProjMPS, y::ProjMPS; kwargs...) = Base.isapprox(x.data, y.data, kwargs...)

function isprojectedat(obj::ProjMPS, ind::IndsT)::Bool where {IndsT}
    return isprojectedat(obj.projector, ind)
end

function Base.:+(a::ProjMPS, b::ProjMPS)::ProjMPS
    return ProjMPS(_add_directsum(a.data, b.data), a.projector | b.projector)
end

function Base.:*(a::ProjMPS, b::Number)::ProjMPS
    return ProjMPS(a.data * b, a.projector)
end

function Base.:*(a::Number, b::ProjMPS)::ProjMPS
    return ProjMPS(b.data * a, b.projector)
end

function Base.:-(obj::ProjMPS)::ProjMPS
    return ProjMPS(-1 * obj.data, obj.projector)
end

function ITensors.truncate(obj::ProjMPS; kwargs...)::ProjMPS
    return project(ProjMPS(truncate(obj.data; kwargs...)), obj.projector)
end
