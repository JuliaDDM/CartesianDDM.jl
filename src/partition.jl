"""
    CartesianBooleanPartition

Boolean partition of unity on a Cartesian domain.

"""
struct CartesianBooleanPartition{T,R} <: DDMMatrix{T}
    decomp::R
    part::R
end

CartesianBooleanPartition{T}(decomp::R, part::R) where {T,R} =
    CartesianBooleanPartition{T,R}(decomp, part)

function size(D::CartesianBooleanPartition)
    (; decomp, part) = D
    ntuple(i -> prod(length.(decomp)), 2)
end

function getindex(A::CartesianBooleanPartition{T}, i, j) where {T,P}
    @boundscheck checkbounds(A, i, j)
    (; decomp, part) = A
    i ≠ j && return zero(T)
    index = getindex(CartesianIndices(decomp), i)
    in(index, CartesianIndices(part)) ? one(T) : zero(T)
end

