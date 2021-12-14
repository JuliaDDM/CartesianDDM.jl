const DDMMatrix{T} = DDMArray{T,2}

"""
    CartesianDDMMatrix{T,C,A}

Objects of `CartesianDDMMatrix` behave like a normal matrix, except
that the storage is not stored continuously but by DDM subdomain.

For now, only Dirichlet matrices are considered.

"""
struct CartesianDDMMatrix{T,C,A} <: DDMMatrix{T}
    context::C
    parent::A
end

function size(x::CartesianDDMMatrix)
    (; context, parent) = x
    (; dims) = context
    ntuple(i -> mapreduce(getdof, *, dims), 2)
end

"""
    getindex(x::CartesianDDMMatrix, i)

Decomposed matrix along rows.
Block-diagonal matrix.

"""
function getindex(x::CartesianDDMMatrix, i::Int, j::Int)
    @boundscheck checkbounds(x, i, j)

    (; context, parent) = x

    (; dims) = context
    gi = getindex(CartesianIndices(getdof.(dims)), i)

    nocontext = removeoverlap(context)
    noglb, _ = ranges(nocontext)

    index = findfirst(noglb) do el
        in(gi, CartesianIndices(el))
    end

    glb, _, lcl = ranges(context)
    domain = CartesianIndices(glb[index])

    li = findfirst(domain) do li
        li == gi
    end

    ii = LinearIndices(CartesianIndices(lcl[index]))[li]

    gj = getindex(CartesianIndices(getdof.(dims)), j)

    lj = findfirst(domain) do lj
        lj == gj
    end

    # something fishy when over = 1
    isnothing(lj) && return zero(eltype(x))

    jj = LinearIndices(CartesianIndices(lcl[index]))[lj]

    parent[index][ii, jj]
end

function decompose(context::CartesianDDMContext, a::AbstractMatrix)
    (; dims) = context

    glb, _ = ranges(context)

    cindices = CartesianIndices(getdof.(dims))
    lindices = LinearIndices(cindices)

    parent = map(glb) do el
        domain = reshape(lindices[CartesianIndices(el)], :)
        a[domain, domain]
    end

    CartesianDDMMatrix{eltype(a),typeof(context),typeof(parent)}(context, parent)
end

import Base.*

function (*)(A::CartesianDDMMatrix, x::CartesianDDMVector)
    context = A.context
    @assert context == x.context

    parent = A.parent .* x.parent

    T = eltype(eltype(parent))

    CartesianDDMVector{T,typeof(context),typeof(parent)}(context, parent)
end

