"""
Add multiple ProjMPS objects on the same projector.

If the bond dimension of the result reaches `maxdim`,
perform patching recursively to reduce the bond dimension.
"""
function _add_patching(
    prjmpss::AbstractVector{ProjMPS};
    cutoff=0.0,
    maxdim=typemax(Int),
    alg=ITensors.Algorithm"fit"(),
    patchorder=Index[],
)::Vector{ProjMPS}
    if length(unique([prjmps.projector for prjmps in prjmpss])) != 1
        error("All ProjMPS objects must have the same projector.")
    end

    # First perform addition upto given maxdim
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
