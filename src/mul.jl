"""
By default, elementwise multiplication will be performed.
"""
function Quantics.automul(
    M1::BlockedMPS,
    M2::BlockedMPS;
    tag_row::String="",
    tag_shared::String="",
    tag_col::String="",
    alg="naive",
    kwargs...,
)
    length.(siteinds(M1)) .== 1 || error("M1 should have only 1 site index per site")
    length.(siteinds(M2)) .== 1 || error("M2 should have only 1 site index per site")
    
    sites_row = _findallsiteinds_by_tag(M1; tag=tag_row)
    sites_shared = _findallsiteinds_by_tag(M1; tag=tag_shared)
    sites_col = _findallsiteinds_by_tag(M2; tag=tag_col)
    sites_matmul = Set(Iterators.flatten([sites_row, sites_shared, sites_col]))

    sites1_ewmul = setdiff(only.(siteinds(M1)), sites_matmul)
    sites2_ewmul = setdiff(only.(siteinds(M2)), sites_matmul)
    sites2_ewmul == sites1_ewmul || error("Invalid sites for elementwise multiplication")

    M1 = Quantics.makesitediagonal(M1, sites1_ewmul; plev=1)
    M2 = Quantics.makesitediagonal(M2, sites2_ewmul; plev=0)

    M1 = Quantics.rearrange_siteinds(
        M1, combinesites(deepcopy(siteinds(M1)), sites_row, sites_shared)
    )

    M2 = Quantics.rearrange_siteinds(
        M2, combinesites(deepcopy(siteinds(M2)), sites_shared, sites_col)
    )

    M = contract(M1, M2; alg=alg, kwargs...)

    M = Quantics.extractdiagonal(M, sites1_ewmul)

    return Quantics.rearrange_siteinds(
        M, collect(Iterators.flatten(siteinds.(M)))
    )
end


function combinesites(sites::Vector{Vector{Index}}, site1::AbstractVector{Index}, site2::AbstractVector{Index})
    length(site1) == length(site2) || error("Length mismatch")
    for (s1, s2) in zip(site1, site2)
        sites = combinesites(sites, s1, s2)
    end
    return sites
end

function combinesites(sites::Vector{Vector{Index}}, site1::Index, site2::Index)
    sites = deepcopy(sites)
    p1 = findfirst(x=>x[1] == site1, sites)
    p2 = findfirst(x=>x[1] == site2, sites)
    if p1 === nothing || p2 === nothing
        error("Site not found")
    end
    if abs(p1 - p2) != 1
        error("Sites are not adjacent")
    end
    deleteat!(sites, minimum(p1, p2))
    deleteat!(sites, minimum(p1, p2))
    insert!(sites, minimum(p1, p2), [site1, site2])
    return sites
end

_findallsiteinds_by_tag(M::BlockMPS; tag=tag) = Quantics.findallsiteinds_by_tag(only.(siteinds(M)); tag=tag)

#==
function combinesites(M::ProjMPS, site1::Index, site2::Index)::ProjMPS
    return project(Quantics.combinesites(Quantics.asMPO(MPS(M)), site1, site2), M.projector)
end

function combinesites(M::BlockedMPS, site1::Index, site2::Index)::BlockedMPS
    return BlockedMPS([combinesites(projmps, site1, site2) for projmps in M])
end

function preprocess(mul::Quantics.MatrixMultiplier{T}, M1::BlockedMPS, M2::BlockedMPS) where {T}
    for (site_row, site_shared, site_col) in
        zip(mul.sites_row, mul.sites_shared, mul.sites_col)
        M1, M2 = combinesites(M1, site_row, site_shared),
        combinesites(M2, site_col, site_shared)
    end
    return M1, M2
end

function postprocess(mul::Quantics.MatrixMultiplier{T}, M::BlockedMPS)::BlockedMPS where {T}
    return BlockedMPS([Quantics.postprocess(mul, x) for x in M])
end

function _postprocess(mul::Quantics.MatrixMultiplier{T}, M::ProjMPS)::ProjMPS where {T}
    return Quantics.postprocess(mul, Quanitcs.asMPO(MPS(M)) )
end

function _preprocess1(mul::Quantics.ElementwiseMultiplier{T}, M1::ProjMPS) where {T}
    return project(_preprocess1(mul, Quantics.asMPO(M1.data)), M1.projector)
end

function _preprocess2(mul::Quantics.ElementwiseMultiplier{T}, M2::ProjMPS) where {T}
    return project(_preprocess2(mul, Quantics.asMPO(M2.data)), M2.projector)
end

function _preprocess1(mul::Quantics.ElementwiseMultiplier{T}, M1::MPO) where {T}
    tensors1 = ITensors.data(M1)
    for s in mul.sites
        p = findfirst(hasind(s), tensors1)
        tensors1[p] = Quantics._asdiagonal(tensors1[p], s)
        replaceind!(tensors1[p], s' => s'')
        replaceind!(tensors1[p], s => s')
    end
    return MPO(tensors1)
end

function _preprocess2(mul::Quantics.ElementwiseMultiplier{T}, M2::MPO) where {T}
    tensors2 = ITensors.data(M2)
    for s in mul.sites
        p = findfirst(hasind(s), tensors2)
        tensors2[p] = Quantics._asdiagonal(tensors2[p], s)
    end
    return MPO(tensors2)
end

function preprocess(mul::Quantics.ElementwiseMultiplier{T}, M1::BlockedMPS, M2::BlockedMPS) where {T}
    return BlockedMPS([_preprocess1(mul, x) for x in M1]), BlockedMPS([_preprocess2(mul, x) for x in M2])
end

function postprocess(mul::Quantics.ElementwiseMultiplier{T}, M::BlockedMPS)::BlockedMPS where {T}
    return BlockedMPS([Quantics.postprocess(mul, x) for x in M])
end
==#