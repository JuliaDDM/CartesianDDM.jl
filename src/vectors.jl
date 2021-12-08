const DDMVector{T} = DDMArray{T,1}

"""
    CartesianDDMVector{T,D,P,A}

# Note on `DenseArray`s

1. Would `CartesianDDMVector` qualify as a `DenseVector` if storage were contiguous?
1. Optional field to signal/ensure coherence?

"""
struct CartesianDDMVector{T,P,A} <: DDMVector{T}
    pou::P
    parent::A
end

function CartesianDDMVector{T}(::UndefInitializer, decomp, part) where {T}
    parent = map(decomp) do indices
        Vector{T}(undef, prod(length.(indices)))
    end

    pou = map(decomp, part) do indices...
        CartesianBooleanPartition{T}(indices...)
    end

    CartesianDDMVector{T,typeof(pou),typeof(parent)}(pou, parent)
end

function CartesianDDMVector(init::Function, decomp, part)
    parent = map(decomp) do indices
        init(prod(length.(indices)))
    end

    T = eltype(eltype(parent))

    pou = map(decomp, part) do indices...
        CartesianBooleanPartition{T}(indices...)
    end

    CartesianDDMVector{T,typeof(pou),typeof(parent)}(pou, parent)
end

function size(x::CartesianDDMVector)
    (; pou, parent) = x
    mapreduce(+, pou) do D
        (; decomp, part) = D
        prod(length.(part))
    end |> tuple
end

"""
    getindex(x::CartesianDDMVector, i)

Assumes partition of unity is Boolean for now. This is a very inefficient
way to access data. It is included only for debugging purposes.

"""
function getindex(x::CartesianDDMVector, i)
    @boundscheck checkbounds(x, i)

    (; pou, parent) = x

    parts = map(pou) do D
        getproperty(D, :part)
    end

    decomps = map(pou) do D
        getproperty(D, :decomp)
    end

    indices = map(first(parts), last(parts)) do start, stop
        first(start):last(stop)
    end

    ci = CartesianIndices(indices)[i]

    index = findfirst(parts) do part
        in(ci, CartesianIndices(part))
    end

    cj = findfirst(CartesianIndices(decomps[index])) do cj
        ci == cj
    end

    j = LinearIndices(CartesianIndices(decomps[index]))[cj]

    parent[index][j]
end

"""
    setindex!(x::CartesianDDMVector, val, i)

This is a very inefficient way to access data. It is included only
for debugging purposes.

"""
function setindex!(x::CartesianDDMVector, val, i)
    @boundscheck checkbounds(x, i)

    (; pou, parent) = x

    parts = map(pou) do D
        getproperty(D, :part)
    end

    decomps = map(pou) do D
        getproperty(D, :decomp)
    end

    indices = map(first(parts), last(parts)) do start, stop
        first(start):last(stop)
    end

    ci = CartesianIndices(indices)[i]

    index = findfirst(parts) do part
        in(ci, CartesianIndices(part))
    end

    cj = findfirst(CartesianIndices(decomps[index])) do cj
        ci == cj
    end

    j = LinearIndices(CartesianIndices(decomps[index]))[cj]

    parent[index][j] = val
end

#=
function CartesianDDMVector{T}(::UndefInitializer, decomp) where {T}
    parent = map(decomp) do indices
        Vector{T}(undef, prod(length.(indices)))
    end
    CartesianDDMVector{T,typeof(parent)}(parent)
end
=#
#=
function CartesianDDMVector(init::Function, decomp)
    parent = map(decomp) do indices
        init(prod(length.(indices)))
    end
    CartesianDDMVector{eltype(first(parent)),typeof(parent)}(parent)
end

=#

