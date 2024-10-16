"""
An MPS with a projector.
"""
struct ProjMPS
    data::MPS
    projector::Projector

    function ProjMPS(data::AbstractMPS, projector::Projector)
        _check_projector_compatibility(projector, data) || error(
            "Incompatible projector and data. Even small numerical noise can cause this error.",
        )
        return new(MPS([x for x in data]), projector)
    end
end

ITensors.siteinds(obj::ProjMPS) = siteinds(obj.data)

function ProjMPS(Ψ::AbstractMPS)
    IndsT = eltype(siteinds(Ψ))
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
        MPS([project(projΨ.data[n], newprj) for n in 1:length(projΨ.data)]), newproj
    )
end

function _iscompatible(projector::Projector, tensor::ITensor)
    # Lazy impmentation
    return norm(project(tensor, projector) - tensor) == 0.0
end

function _check_projector_compatibility(
    projector::Projector, Ψ::MPS, sites::AbstractVector{<:AbstractVector}
)
    links = linkinds(Ψ)
    sitedims = [collect(dim.(s)) for s in sites]

    sitetensors = []
    push!(
        sitetensors,
        reshape(
            Array(Ψ[1], [sites[1]..., links[1]]), [1, prod(sitedims[1]), dim(links[1])]...
        ),
    )
    for n in 2:(length(Ψ) - 1)
        push!(
            sitetensors,
            reshape(
                Array(Ψ[n], [links[n - 1], sites[n]..., links[n]]),
                dim(links[n - 1]),
                prod(sitedims[n]),
                dim(links[n]),
            ),
        )
    end
    push!(
        sitetensors,
        reshape(
            Array(Ψ[end], [links[end], sites[end]...]),
            dim(links[end]),
            prod(sitedims[end]),
            1,
        ),
    )

    return reduce(
        &,
        _check_projector_compatibility(projector[n], sitedims[n], sitetensors[n]) for
        n in 1:length(Ψ)
    )
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
    for (p, s) in zip(projmps.projector, projmps.sites)
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
    sites_diagonal = siteinds(all, mps_diagonal)
    projmps_diagonal = ProjMPS(mps_diagonal, sites_diagonal)

    target_positions = Quantics.findallsiteinds_by_tag(siteinds(MPS(projmps)); tag=tag)
    prjsiteinds = Dict{Index{Int},Int}()
    for (p, s) in zip(projmps.projector, projmps.sites)
        for (p_, s_) in zip(p, s)
            iszero(p_) && continue
            prjsiteinds[s_] = p_
            if s_ ∈ target_positions
                prjsiteinds[s_'] = p_
            end
        end
    end

    return project(projmps_diagonal, prjsiteinds)
end

function Quantics.makesitediagonal(projmpss::ProjMPSContainer, sites)
    return ProjMPSContainer([
        Quantics.makesitediagonal(projmps, sites) for projmps in projmpss.data
    ])
end

function Quantics.extractdiagonal(projmps::ProjMPS, tag::String)
    mps_diagonal = Quantics.extractdiagonal(MPS(projmps), tag)
    sites_diagonal = siteinds(all, mps_diagonal)
    projmps_diagonal = ProjMPS(mps_diagonal, sites_diagonal)
    sites_diagonal_set = Set(Iterators.flatten(sites_diagonal))

    prjsiteinds = Dict{Index{Int},Int}()
    for (p, s) in zip(projmps.projector, projmps.sites)
        for (p_, s_) in zip(p, s)
            !iszero(p_) || continue
            s_ ∈ sites_diagonal_set || continue
            prjsiteinds[s_] = p_
        end
    end

    return project(projmps_diagonal, prjsiteinds)
end

function Quantics.extractdiagonal(projmpss::ProjMPSContainer, sites)
    return ProjMPSContainer([
        Quantics.extractdiagonal(projmps, sites) for projmps in projmpss.data
    ])
end

function Quantics.rearrange_siteinds(projmps::ProjMPS, sites)
    mps_rearranged = Quantics.rearrange_siteinds(MPS(projmps), sites)
    projmps_rearranged = ProjMPS(mps_rearranged, sites)
    prjsiteinds = Dict{Index{Int},Int}()
    for (p, s) in zip(projmps.projector, projmps.sites)
        for (p_, s_) in zip(p, s)
            if p_ != 0
                prjsiteinds[s_] = p_
            end
        end
    end
    return project(projmps_rearranged, prjsiteinds)
end

function Quantics.rearrange_siteinds(projmpss::ProjMPSContainer, sites)
    return ProjMPSContainer([
        Quantics.rearrange_siteinds(projmps, sites) for projmps in projmpss.data
    ])
end

# Miscellaneous Functions
function Base.show(io::IO, obj::ProjMPS)
    return print(io, "ProjMPS projected on $(obj.projector.data)")
end

function ITensors.prime(Ψ::ProjMPS, args...; kwargs...)
    return ProjMPS(
        prime(MPS(Ψ), args...; kwargs...), prime.(Ψ.sites, args...; kwargs...), Ψ.projector
    )
end

function ITensors.prime(Ψ::ProjMPSContainer, args...; kwargs...)
    return ProjMPSContainer([prime(projmps, args...; kwargs...) for projmps in Ψ.data])
end

Base.isapprox(x::ProjMPS, y::ProjMPS; kwargs...) = Base.isapprox(x.data, y.data, kwargs...)

# Random MPO Functions (commented out)
#==
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
    return MPS(tensors)
end

function _random_mpo(sites::AbstractVector{<:AbstractVector{Index{T}}}; m::Int=1) where {T}
    return _random_mpo(Random.default_rng(), sites; m=m)
end
==#
