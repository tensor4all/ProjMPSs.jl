function ITensors.contract(M1::ProjMPS, M2::ProjMPS; kwargs...)::Union{ProjMPS,Nothing}
    if !hasoverlap(M1.projector, M2.projector)
        return nothing
    end
    proj, _ = _projector_after_contract(M1, M2)

    Ψ = FMPOC.contract_mpo_mpo(MPO(collect(M1.data)), MPO(collect(M2.data)); kwargs...)
    return project(ProjMPS(Ψ), proj)
end

# Figure out `projector` after contracting ProjMPS objects
function _projector_after_contract(M1::ProjMPS, M2::ProjMPS)
    sites1 = _allsites(M1)
    sites2 = _allsites(M2)

    external_sites = setdiff(union(sites1, sites2), intersect(sites1, sites2))

    proj = deepcopy(M1.projector.data)
    empty!(proj)

    for s in external_sites
        if isprojectedat(M1, s)
            proj[s] = M1.projector[s]
        end
        if isprojectedat(M2, s)
            proj[s] = M2.projector[s]
        end
    end

    return Projector(proj), external_sites
end

function _is_externalsites_compatible_with_projector(external_sites, projector)
    for s in keys(projector)
        if !(s ∈ external_sites)
            return false
        end
    end
    return true
end

"""
Project two ProjMPS objects to `proj` before contracting them.
"""
function projcontract(
    M1::ProjMPS, M2::ProjMPS, proj::Projector; kwargs...
)::Union{Nothing,ProjMPS}
    # Project M1 and M2 to `proj` before contracting
    M1 = project(M1, proj)
    M2 = project(M2, proj)
    if M1 === nothing || M2 === nothing
        return nothing
    end

    _, external_sites = _projector_after_contract(M1, M2)

    if !_is_externalsites_compatible_with_projector(external_sites, proj)
        error("The projector contains projection onto a site what is not a external sites.")
    end

    return ITensors.contract(M1, M2; kwargs...)
end

"""
Project two ProjMPS objects to `proj` before contracting them.
The results are summed.
"""
function projcontract(
    M1::AbstractVector{ProjMPS},
    M2::AbstractVector{ProjMPS},
    proj::Projector;
    alg::String="naive",
    kwargs...,
)::Union{Nothing,ProjMPS}
    results = ProjMPS[]
    for M1_ in M1
        for M2_ in M2
            r = projcontract(M1_, M2_, proj; alg, kwargs...)
            if r !== nothing
                push!(results, r)
            end
        end
    end

    if isempty(results)
        return nothing
    end

    if length(results) == 1
        return results[1]
    end

    res = _add(results; kwargs...)

    return res
end

_add_directsum(a::AbstractMPS, b::AbstractMPS) = +(a, b; alg="directsum")

function _add(Ψs::AbstractVector{MPO}; cutoff=0.0, maxdim=typemax(Int), nsweeps=2)::MPO
    Ψsum::MPO = deepcopy(Ψs[1])
    for Ψ::MPO in Ψs[2:end]
        Ψsum = +(Ψsum, Ψ; alg="directsum")
        truncate!(Ψsum; cutoff, maxdim)
    end
    return FMPOC.fit(Ψs, Ψsum; maxdim, cutoff, nsweeps)
end

function _add(
    Ψs::AbstractVector{ProjMPS}; cutoff=0.0, maxdim=typemax(Int), nsweeps=2
)::ProjMPS
    mpos::Vector{MPO} = [MPO(collect(x.data)) for x in Ψs]
    sum_mps = _add(mpos; cutoff, maxdim, nsweeps)
    newprj = reduce(|, (x.projector for x in Ψs))
    return project(ProjMPS(sum_mps), newprj)
end

function ITensors.contract(M1::BlockedMPS, M2::BlockedMPS; kwargs...)::Union{BlockedMPS}
    blocks = OrderedSet((
        _projector_after_contract(b1, b2)[1] for b1 in M1.data, b2 in M2.data
    ))
    prjmpss = ProjMPS[]
    M1_::Vector{ProjMPS} = collect(M1)
    M2_::Vector{ProjMPS} = collect(M2)
    for b in blocks
        res = projcontract(M1_, M2_, b; kwargs...)
        if res !== nothing
            push!(prjmpss, res)
        end
    end
    return BlockedMPS(prjmpss)
end
