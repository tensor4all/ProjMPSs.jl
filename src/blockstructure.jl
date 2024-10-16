"""
Consists of a collection of non-overlapping blocks.
The blocks do not nessarily have to cover the whole space.

FIXME: We do not need this struc? We can just use a vector of projectors.
"""
struct BlockStructure{IndsT}
    blocks::OrderedSet{Projector{IndsT}} # Non-overlapping blocks

    function BlockStructure(blocks::AbstractVector{Projector{IndsT}}) where {IndsT}
        for (ib1, bl1) in enumerate(blocks), (ib2, bl2) in enumerate(blocks)
            if ib1 != ib2
                if hasoverlap(bl1, bl2)
                    error("Blocks are overlapping: $bl1 and $bl2")
                end
            end
        end
        return new{IndsT}(OrderedSet(blocks))
    end
end

Base.length(bs::BlockStructure) = length(bs.blocks)

function Base.iterate(p::BlockStructure, state=1)
    if state > length(p.blocks)
        return nothing
    end
    return (p.blocks[state], state + 1)
end

Base.getindex(p::BlockStructure, index::Int) = p.blocks[index]
