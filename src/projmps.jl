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

function find_nested_index(data::Vector{Vector{T}}, target::T) where {T}
    for (i, subvector) in enumerate(data)
        j = findfirst(x -> x == target, subvector)
        if j !== nothing
            return (i, j)
        end
    end
    return nothing  # Not found
end

# Quantics Functions
function Quantics.makesitediagonal(projmps::ProjMPS, site::Index)
    mps_diagonal = Quantics.makesitediagonal(MPS(projmps), site)
    sites_diagonal = siteinds(all, mps_diagonal)
    projmps_diagonal = ProjMPS(mps_diagonal, sites_diagonal)

    prjsiteinds = Dict{Index{Int},Int}()
    for (p, s) in zip(projmps.projector, siteinds(projmps))
        for (p_, s_) in zip(p, s)
            iszero(p_) && continue
            prjsiteinds[s_] = p_
            if s_ == site
                prjsiteinds[s_'] = p_
            end
        end
    end

    return project(projmps_diagonal, prjsiteinds)
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

function Quantics.extractdiagonal(projmps::ProjMPS, tag::String)
    mps_diagonal = Quantics.extractdiagonal(MPS(projmps), tag)
    projmps_diagonal = ProjMPS(mps_diagonal)

    target_sites = Quantics.findallsiteinds_by_tag(
        unique(noprime.(Iterators.flatten(siteinds(projmps)))); tag=tag
    )

    newproj = deepcopy(projmps.projector)

    for s in target_sites
        if isprojectedat(newproj, s)
            delete!(newproj, s')
        end
    end

    return project(projmps_diagonal, newproj)
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
