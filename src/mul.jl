"""
By default, elementwise multiplication will be performed.
"""
function automul(
    M1::BlockedMPS,
    M2::BlockedMPS;
    tag_row::String="",
    tag_shared::String="",
    tag_col::String="",
    alg="naive",
    kwargs...,
)
    sites_row = Quantics.findallsiteinds_by_tag(siteinds(M1); tag=tag_row)
    sites_shared = Quantics.findallsiteinds_by_tag(siteinds(M1); tag=tag_shared)
    sites_col = Quantics.findallsiteinds_by_tag(siteinds(M2); tag=tag_col)
    sites_matmul = Set(Iterators.flatten([sites_row, sites_shared, sites_col]))

    if sites_shared != Quantics.findallsiteinds_by_tag(siteinds(M2); tag=tag_shared)
        error("Invalid shared sites for MatrixMultiplier")
    end

    matmul = Quantics.MatrixMultiplier(sites_row, sites_shared, sites_col)
    ewmul = Quantics.ElementwiseMultiplier([s for s in siteinds(M1) if s ∉ sites_matmul])

    M1_, M2_ = _preprocess(matmul, M1_, M2_)
    M1_, M2_ = _preprocess(ewmul, M1_, M2_)

    M = contract(M1_, M2_; alg=alg, kwargs...)

    M = _postprocess(matmul, M)
    M = _postprocess(ewmul, M)

    if in(:maxdim, keys(kwargs))
        return truncate(M; maxdim=kwargs[:maxdim])
    else
        return M
    end
end

function combinesites(M::BlockedMPS, site1::Index, site2::Index)::BlockedMPS
    return BlockedMPS([combine_sites(projmps, site1, site2) for projmps in M])
end

function preprocess(mul::MatrixMultiplier{T}, M1::BlockedMPS, M2::BlockedMPS) where {T}
    for (site_row, site_shared, site_col) in
        zip(mul.sites_row, mul.sites_shared, mul.sites_col)
        M1, M2 = combinesites(M1, site_row, site_shared),
        combinesites(M2, site_col, site_shared)
    end
    return M1, M2
end

function postprocess(mul::MatrixMultiplier{T}, M::BlockedMPS)::BlockedMPS where {T}
    return BlockedMPS([Quantics.postprocess(mul, x) for x in M])
end

function _preprocess1(mul::ElementwiseMultiplier{T}, M1::MPO) where {T}
    tensors1 = ITensors.data(M1)
    for s in mul.sites
        p = findfirst(hasind(s), tensors1)
        tensors1[p] = _asdiagonal(tensors1[p], s)
        replaceind!(tensors1[p], s' => s'')
        replaceind!(tensors1[p], s => s')
    end
    return MPO(tensors1)
end

function _preprocess2(mul::ElementwiseMultiplier{T}, M2::MPO) where {T}
    tensors2 = ITensors.data(M2)
    for s in mul.sites
        p = findfirst(hasind(s), tensors2)
        tensors2[p] = _asdiagonal(tensors2[p], s)
    end
    return MPO(tensors2)
end

function preprocess(mul::ElementwiseMultiplier{T}, M1::BlockedMPS, M2::BlockedMPS) where {T}
    return [_preprocess1(mul, x) for x in M1], [_preprocess2(mul, x) for x in M2]
end

function postprocess(mul::ElementwiseMultiplier{T}, M::BlockedMPS)::BlockedMPS where {T}
    return BlockedMPS([Quantics.postprocess(mul, x) for x in M])
end