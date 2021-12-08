"""
    CartesianBooleanPartition

Boolean partition of unity on a Cartesian domain.

"""
struct CartesianBooleanPartition{T,R} <: DDMMatrix{T}
    part::R
end

CartesianBooleanPartition{T}(part) where {T} =
    CartesianBooleanPartition{T,typeof(part)}(part)

