"""
    CartesianBooleanPartition

Boolean partition of unity on a Cartesian domain.

1. Can `eltype` be `Bool`?

"""
struct CartesianBooleanPartition{T,R} <: DDMMatrix{T}
    decomp::R
    part::R
end

CartesianBooleanPartition{T}(decomp::R, part::R) where {T,R} =
    CartesianBooleanPartition{T,R}(decomp, part)

function size(D::CartesianBooleanPartition)
    (; decomp, part) = D

    n = mapreduce(+, decomp) do indices
        prod(length.(indices))
    end

    ntuple(i -> n, 2)
end

function getindex(A::CartesianBooleanPartition{T}, i::Int, j::Int) where {T}
    @boundscheck checkbounds(A, i, j)

    (; decomp, part) = A

    i != j && return zero(T)

    (; indices, nproc, nover) = decomp

    subindices = OneTo.(nover .* (nproc .- 1) .+ length.(indices))

    ci = CartesianIndices(subindices)[i]

    @show ci

#=
    index = findfirst(decomp) do sub
        in(ci, CartesianIndices(sub))
    end

    if in(ci, CartesianIndices(part[index]))
        return one(T)
    else
        return zero(T)
    end
    =#
    zero(T)
end

