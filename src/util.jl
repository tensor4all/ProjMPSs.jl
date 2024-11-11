
function allequal(collection)
    if isempty(collection)
        return true
    end
    c = first(collection)
    return all([x == c for x in collection])
end

function Not(index::Int, length::Int)
    return vcat(1:(index - 1), (index + 1):length)
end

function typesafe_iterators_product(::Val{N}, dims) where {N}
    return Iterators.product(ntuple(i -> 1:dims[i], N)...)
end

_getindex(x, indices) = ntuple(i -> x[indices[i]], length(indices))

function _contract(
    a::AbstractArray{T1,N1},
    b::AbstractArray{T2,N2},
    idx_a::NTuple{n1,Int},
    idx_b::NTuple{n2,Int},
) where {T1,T2,N1,N2,n1,n2}
    length(idx_a) == length(idx_b) || error("length(idx_a) != length(idx_b)")
    # check if idx_a contains only unique elements
    length(unique(idx_a)) == length(idx_a) || error("idx_a contains duplicate elements")
    # check if idx_b contains only unique elements
    length(unique(idx_b)) == length(idx_b) || error("idx_b contains duplicate elements")
    # check if idx_a and idx_b are subsets of 1:N1 and 1:N2
    all(1 <= idx <= N1 for idx in idx_a) || error("idx_a contains elements out of range")
    all(1 <= idx <= N2 for idx in idx_b) || error("idx_b contains elements out of range")

    rest_idx_a = setdiff(1:N1, idx_a)
    rest_idx_b = setdiff(1:N2, idx_b)

    amat = reshape(
        permutedims(a, (rest_idx_a..., idx_a...)),
        prod(_getindex(size(a), rest_idx_a)),
        prod(_getindex(size(a), idx_a)),
    )
    bmat = reshape(
        permutedims(b, (idx_b..., rest_idx_b...)),
        prod(_getindex(size(b), idx_b)),
        prod(_getindex(size(b), rest_idx_b)),
    )

    return reshape(
        amat * bmat, _getindex(size(a), rest_idx_a)..., _getindex(size(b), rest_idx_b)...
    )
end

# QUESTION: Is this really a shallow copy? It works like a deep copy - Gianluca 
function shallowcopy(original)
    fieldnames = Base.fieldnames(typeof(original))
    new_fields = [Base.copy(getfield(original, name)) for name in fieldnames]
    return (typeof(original))(new_fields...)
end

function _asdiagonal(t, site::Index{T}; baseplev=0)::ITensor where {T<:Number}
    ITensors.hasinds(t, site') && error("Found $(site')")
    links = ITensors.uniqueinds(ITensors.inds(t), site)
    rawdata = Array(t, links..., site)
    tensor = zeros(eltype(t), size(rawdata)..., ITensors.dim(site))
    for i in 1:ITensors.dim(site)
        tensor[.., i, i] = rawdata[.., i]
    end
    return ITensor(
        tensor, links..., ITensors.prime(site, baseplev + 1), ITensors.prime(site, baseplev)
    )
end
