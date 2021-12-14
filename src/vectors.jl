const DDMVector{T} = DDMArray{T,1}

"""
    CartesianDDMVector{T,C,A}

Objects of `CartesianDDMVector` behave like a normal vector, except
that the storage is not stored continuously but by DDM subdomain.

For now, a Boolean partition of unity is assumed but this could be
generalized at a later time.

# Note on `DenseArray`s

1. Would `CartesianDDMVector <: DenseVector` if contiguous storage?
1. Optional field to signal/ensure coherence?

"""
struct CartesianDDMVector{T,C,A} <: DDMVector{T}
    context::C
    parent::A
end

function size(x::CartesianDDMVector)
    (; context, parent) = x
    (; dims) = context
    (mapreduce(getdof, *, dims),)
end

"""
    getindex(x::CartesianDDMVector, i)

Assumes partition of unity is Boolean for now. This is a very
inefficient way to read data. It is included only for debugging
purposes.

"""
function getindex(x::CartesianDDMVector, i::Int)
    @boundscheck checkbounds(x, i)

    (; context, parent) = x

    (; dims) = context
    ci = getindex(CartesianIndices(getdof.(dims)), i)

    nocontext = removeoverlap(context)
    noglb, _ = ranges(nocontext)

    index = findfirst(noglb) do el
        in(ci, CartesianIndices(el))
    end

    glb, _, lcl = ranges(context)
    domain = CartesianIndices(glb[index])

    cj = findfirst(domain) do cj
        cj == ci
    end

    j = LinearIndices(CartesianIndices(lcl[index]))[cj]

    parent[index][j]
end

"""
    setindex!(x::CartesianDDMVector, val, i::Int)

Writes data in every subdomain that contains element i. This is
a very inefficient way to write data. It is included only for
debugging purposes.

"""
function setindex!(x::CartesianDDMVector, val, i::Int)
    @boundscheck checkbounds(x, i)

    (; context, parent) = x

    (; dims) = context
    ci = getindex(CartesianIndices(getdof.(dims)), i)

    glb, _, lcl = ranges(context)

    indices = findall(glb) do el
        in(ci, CartesianIndices(el))
    end

    for index in indices
        domain = CartesianIndices(glb[index])

        cj = findfirst(domain) do cj
            cj == ci
        end

        domain = CartesianIndices(lcl[index])

        j = LinearIndices(domain)[cj]

        parent[index][j] = val
    end

    val
end

#=
"""
Do not use until `makecoherent` implemented.

"""
function CartesianDDMVector{T}(::UndefInitializer, context) where {T}
    _, lcl = ranges(context)

    parent = map(lcl) do el
        Array{T}(undef, prod(length.(el)))
    end

    CartesianDDMVector{T,typeof(context),typeof(parent)}(context, parent)
end

"""
Do not use until `makecoherent` implemented.

"""
function CartesianDDMVector(init::Function, context)
    _, lcl = ranges(context)

    parent = map(lcl) do el
        init(prod(length.(el)))
    end

    T = eltype(eltype(parent))

    CartesianDDMVector{T,typeof(context),typeof(parent)}(context, parent)
end
=#

function decompose(context::CartesianDDMContext, x::AbstractVector)
    (; dims) = context

    y = reshape(x, getdof.(dims))

    glb, _ = ranges(context)

    parent = map(glb) do el
        y[CartesianIndices(el)]
    end

    CartesianDDMVector{eltype(x),typeof(context),typeof(parent)}(context, parent)
end

