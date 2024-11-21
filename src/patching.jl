"""
Add multiple ProjMPS objects on the same projector.

If the bond dimension of the result reaches `maxdim`,
perform patching recursively to reduce the bond dimension.
"""
function _add_patching(
    prjmpss::AbstractVector{ProjMPS};
    cutoff=0.0,
    maxdim=typemax(Int),
    alg="fit",
    patchorder=Index[],
)::Vector{ProjMPS}
    if length(unique([prjmps.projector for prjmps in prjmpss])) != 1
        error("All ProjMPS objects must have the same projector.")
    end

    # First perform addition upto given maxdim
    # TODO: Early termination if the bond dimension reaches maxdim
    sum_approx = _add(prjmpss...; alg, cutoff, maxdim)

    # If the bond dimension is less than maxdim, return the result
    maxbonddim(sum_approx) < maxdim && return [sum_approx]

    @assert maxbonddim(sum_approx) == maxdim

    nextprjidx = _next_projindex(prjmpss[1].projector, patchorder)

    nextprjidx === nothing && return BlockedMPS(sum_approx)

    blocks = ProjMPS[]
    for prjval in 1:ITensors.dim(nextprjidx)
        prj_ = prjmpss[1].projector & Projector(nextprjidx => prjval)
        blocks =
            blocks ∪ _add_patching(
                [project(prjmps, prj_) for prjmps in prjmpss];
                cutoff,
                maxdim,
                alg,
                patchorder,
            )
    end

    return blocks
end

"""
Return the next index to be projected.
"""
function _next_projindex(prj::Projector, patchorder)::Union{Nothing,Index}
    idx = findfirst(x -> !isprojectedat(prj, x), patchorder)
    if idx === nothing
        return nothing
    else
        return patchorder[idx]
    end
end

"""
Add multiple BlockedMPS objects.
"""
function add_patching(
    bmpss::AbstractVector{BlockedMPS};
    cutoff=0.0,
    maxdim=typemax(Int),
    alg="fit",
    patchorder=Index[],
)::BlockedMPS
    result = _add_patching(union(values(x) for x in bmpss); cutoff, maxdim, alg, patchorder)
    return BlockedMPS(result)
end

"""
Adaptive patching

Do patching recursively to reduce the bond dimension.
If the bond dimension of a ProjMPS exceeds `maxdim`, perform patching.
"""
function adaptive_patching(
    prjmps::ProjMPS, patchorder; cutoff=0.0, maxdim=typemax(Int)
)::Vector{ProjMPS}
    if maxbonddim(prjmps) <= maxdim
        return [prjmps]
    end

    # If the bond dimension exceeds maxdim, perform patching
    refined_prjmpss = ProjMPS[]
    nextprjidx = _next_projindex(prjmps.projector, patchorder)
    if nextprjidx === nothing
        return [prjmps]
    end

    for prjval in 1:ITensors.dim(nextprjidx)
        prj_ = prjmps.projector & Projector(nextprjidx => prjval)
        prjmps_ = truncate(project(prjmps, prj_); cutoff, maxdim)
        if maxbonddim(prjmps_) <= maxdim
            push!(refined_prjmpss, prjmps_)
        else
            append!(refined_prjmpss, adaptive_patching(prjmps_, patchorder; cutoff, maxdim))
        end
    end
    return refined_prjmpss
end

"""
Adaptive patching

Do patching recursively to reduce the bond dimension.
If the bond dimension of a ProjMPS exceeds `maxdim`, perform patching.
"""
function adaptive_patching(
    prjmpss::BlockedMPS, patchorder; cutoff=0.0, maxdim=typemax(Int)
)::BlockedMPS
    return BlockedMPS(
        collect(
            Iterators.flatten((
                apdaptive_patching(prjmps; cutoff, maxdim, patchorder) for
                prjmps in values(prjmpss)
            )),
        ),
    )
end
