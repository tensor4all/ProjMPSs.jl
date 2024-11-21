"""
Conversion from a `TCIAlgorithms.ProjTensorTrain` to a `ProjMPS`.
"""
function ProjMPS(projtt::TCIA.ProjTensorTrain{T}, sites)::ProjMPS where {T}
    links = [Index(ld, "Link,l=$l") for (l, ld) in enumerate(TCI.linkdims(projtt.data))]

    tensors = ITensor[]
    sitedims = [collect(dim.(s)) for s in sites]
    linkdims = dim.(links)

    push!(
        tensors,
        ITensor(
            reshape(projtt.data[1], 1, prod(sitedims[1]), linkdims[1]),
            sites[1]...,
            links[1],
        ),
    )

    for n in 2:(length(projtt.data) - 1)
        push!(
            tensors,
            ITensor(
                reshape(projtt.data[n], linkdims[n - 1], prod(sitedims[n]), linkdims[n]),
                links[n - 1],
                sites[n]...,
                links[n],
            ),
        )
    end

    push!(
        tensors,
        ITensor(
            reshape(projtt.data[end], linkdims[end], prod(sitedims[end])),
            links[end],
            sites[end]...,
        ),
    )

    proj = Dict{Index,Int}()
    for i in eachindex(projtt.projector.data)
        for j in eachindex(projtt.projector.data[i])
            if projtt.projector.data[i][j] > 0
                proj[sites[i][j]] = projtt.projector.data[i][j]
            end
        end
    end

    return ProjMPS(MPS(tensors), Projector(proj))
end

"""
Conversion from a `TCIAlgorithms.Proj` to a `BlockedMPS`.
"""
function BlockedMPS(obj::TCIA.ProjTTContainer{T}, sites)::BlockedMPS where {T}
    return BlockedMPS([ProjMPS(x, sites) for x in obj.data])
end
