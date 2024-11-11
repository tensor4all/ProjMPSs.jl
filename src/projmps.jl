"""
An MPS with a projector.
"""
struct ProjMPS
    data::MPS
    projector::Projector

    function ProjMPS(data::AbstractMPS, projector::Projector)
        _iscompatible(projector, data) || error(
            "Incompatible projector and data. Even small numerical noise can cause this error.",
        )
        projector = _trim_projector(data, projector)
        return new(MPS([x for x in data]), projector)
    end
end

siteinds(obj::ProjMPS) = collect(ITensors.siteinds(MPO([x for x in obj.data])))

_allsites(Ψ::AbstractMPS) = collect(Iterators.flatten(ITensors.siteinds(MPO(collect(Ψ)))))
_allsites(Ψ::ProjMPS) = _allsites(Ψ.data)

maxlinkdim(Ψ::ProjMPS) = ITensors.maxlinkdim(Ψ.data)
maxbonddim(Ψ::ProjMPS) = ITensors.maxlinkdim(Ψ.data)

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
    return ProjMPS(Ψ, Projector())
end

# Conversion Functions
ITensors.MPS(projΨ::ProjMPS) = projΨ.data

function project(tensor::ITensor, projector::Projector)
    slice = Union{Int,Colon}[
        isprojectedat(projector, idx) ? projector[idx] : Colon() for
        idx in ITensors.inds(tensor)
    ]
    data_org = Array(tensor, ITensors.inds(tensor)...)
    data_trim = zero(data_org)
    if all(broadcast(x -> x isa Integer, slice))
        data_trim[slice...] = data_org[slice...]
    else
        data_trim[slice...] .= data_org[slice...]
    end
    return ITensor(data_trim, ITensors.inds(tensor)...)
end

function project(projΨ::ProjMPS, projector::Projector)::Union{Nothing,ProjMPS}
    newprj = projector & projΨ.projector
    if newprj === nothing
        return nothing
    end

    return ProjMPS(
        MPS([project(projΨ.data[n], newprj) for n in 1:length(projΨ.data)]), newprj
    )
end

function project(Ψ::AbstractMPS, projector::Projector)::Union{Nothing,ProjMPS}
    return project(ProjMPS(Ψ), projector)
end

function project(
    projΨ::ProjMPS, projector::Dict{InsT,Int}
)::Union{Nothing,ProjMPS} where {InsT}
    return project(projΨ, Projector(projector))
end

function _iscompatible(projector::Projector, tensor::ITensor)
    # Lazy implementation
    return ITensors.norm(project(tensor, projector) - tensor) == 0.0
end

function _iscompatible(projector::Projector, Ψ::AbstractMPS)
    return all((_iscompatible(projector, x) for x in Ψ))
end

function _makesitediagonal(
    projmps::ProjMPS, sites::AbstractVector{Index{IndsT}}; baseplev=0
) where {IndsT}
    M_ = deepcopy(MPO(collect(MPS(projmps))))
    for site in sites
        target_site::Int = only(ITensors.findsites(M_, site))
        M_[target_site] = _asdiagonal(M_[target_site], site; baseplev=baseplev)
    end
    return project(M_, projmps.projector)
end

function _makesitediagonal(projmps::ProjMPS, site::Index; baseplev=0)
    return _makesitediagonal(projmps, [site]; baseplev=baseplev)
end

function makesitediagonal(projmps::ProjMPS, site::Index)
    return _makesitediagonal(projmps, site; baseplev=0)
end

function makesitediagonal(projmps::ProjMPS, sites::AbstractVector{Index})
    return _makesitediagonal(projmps, sites; baseplev=0)
end

function makesitediagonal(projmps::ProjMPS, tag::String)
    mps_diagonal = Quantics.makesitediagonal(MPS(projmps), tag)
    projmps_diagonal = ProjMPS(mps_diagonal)

    target_sites = Quantics.findallsiteinds_by_tag(
        unique(ITensors.noprime.(Iterators.flatten(siteinds(projmps)))); tag=tag
    )

    newproj = deepcopy(projmps.projector)
    for s in target_sites
        if isprojectedat(projmps.projector, s)
            newproj[ITensors.prime(s)] = newproj[s]
        end
    end

    return project(projmps_diagonal, newproj)
end

# FIXME: may be type unstable
function _find_site_allplevs(tensor::ITensor, site::Index; maxplev=10)
    ITensors.plev(site) == 0 || error("Site index must be unprimed.")
    return [
        ITensors.prime(site, plev) for
        plev in 0:maxplev if ITensors.prime(site, plev) ∈ ITensors.inds(tensor)
    ]
end

function extractdiagonal(
    projmps::ProjMPS, sites::AbstractVector{Index{IndsT}}
) where {IndsT}
    tensors = collect(projmps.data)
    for i in eachindex(tensors)
        for site in intersect(sites, ITensors.inds(tensors[i]))
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

function extractdiagonal(projmps::ProjMPS, tag::String)::ProjMPS
    targetsites = Quantics.findallsiteinds_by_tag(
        unique(ITensors.noprime.(_allsites(projmps))); tag=tag
    )
    return extractdiagonal(projmps, targetsites)
end

function rearrange_siteinds(projmps::ProjMPS, sites)
    mps_rearranged = Quantics.rearrange_siteinds(MPS(projmps), sites)
    return project(ProjMPS(mps_rearranged), projmps.projector)
end

# Miscellaneous Functions
function Base.show(io::IO, obj::ProjMPS)
    return print(io, "ProjMPS projected on $(obj.projector.data)")
end

function prime(Ψ::ProjMPS, args...; kwargs...)
    return ProjMPS(
        ITensors.prime(MPS(Ψ), args...; kwargs...),
        ITensors.prime.(siteinds(Ψ), args...; kwargs...),
        Ψ.projector,
    )
end

Base.isapprox(x::ProjMPS, y::ProjMPS; kwargs...) = Base.isapprox(x.data, y.data, kwargs...)

function isprojectedat(obj::ProjMPS, ind::IndsT)::Bool where {IndsT}
    return isprojectedat(obj.projector, ind)
end

function _fitsum(
    input_states::AbstractVector{T},
    init::T;
    coeffs::AbstractVector{<:Number}=ones(Int, length(input_states)),
    kwargs...,
) where {T}
    if !(:nsweeps ∈ keys(kwargs))
        kwargs = Dict{Symbol,Any}(kwargs)
        kwargs[:nsweeps] = 1
    end
    Ψs = [MPS(collect(x)) for x in input_states]
    init_Ψ = MPS(collect(init))
    res = FMPOC.fit(Ψs, init_Ψ; coeffs=coeffs, kwargs...)
    return T(collect(res))
end

function _add(
    ψ::AbstractMPS...;
    alg=ITensors.Algorithm"fit"(),
    cutoff=1e-15,
    maxdim=typemax(Int),
    kwargs...,
)
    alg = ITensors.Algorithm(alg)
    if alg == ITensors.Algorithm"directsum"()
        return +(ITensors.Algorithm(alg), ψ...)
    elseif alg == ITensors.Algorithm"densitymatrix"()
        if cutoff < 1e-15
            @warn "Cutoff is very small, it may suffer from numerical round errors. The densitymatrix algorithm squares the singular values of the reduce density matrix. Please consider increasing it or using fit algorithm."
        end
        return +(ITensors.Algorithm"densitymatrix"(), ψ...; cutoff, maxdim, kwargs...)
    elseif alg == ITensors.Algorithm"fit"()
        function f(x, y)
            return ITensors.truncate(
                +(ITensors.Algorithm("directsum"), x, y); cutoff, maxdim
            )
        end
        res_dm = reduce(f, ψ)
        res = _fitsum(collect(ψ), res_dm; cutoff, maxdim, kwargs...)
        return res
    else
        error("Unknown algorithm $(alg) for addition!")
    end
end

function Base.:+(
    Ψ::ProjMPS...;
    alg=ITensors.Algorithm"directsum"(),
    cutoff=0.0,
    maxdim=typemax(Int),
    kwargs...,
)::ProjMPS
    return _add(Ψ...; alg=alg, cutoff=cutoff, maxdim=maxdim, kwargs...)
end

function _add(
    Ψ::ProjMPS...;
    alg=ITensors.Algorithm"directsum"(),
    cutoff=0.0,
    maxdim=typemax(Int),
    kwargs...,
)::ProjMPS
    return project(
        _add(
            [x.data for x in Ψ]...;
            alg=ITensors.Algorithm(alg),
            cutoff=cutoff,
            maxdim=maxdim,
        ),
        reduce(|, [x.projector for x in Ψ]),
    )
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

function truncate(obj::ProjMPS; kwargs...)::ProjMPS
    return project(ProjMPS(ITensors.truncate(obj.data; kwargs...)), obj.projector)
end

function _norm(M::AbstractMPS)
    if ITensors.isortho(M)
        return ITensors.norm(M[orthocenter(M)])
    end
    norm2_M = ITensors.dot(M, M)
    return sqrt(abs(norm2_M))
end

function ITensors.norm(M::ProjMPS)
    return _norm(MPS(M))
end
