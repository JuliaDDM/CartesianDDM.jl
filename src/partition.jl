struct BooleanPartition{T,R} <: AbstractMatrix{T}
    part::R
    decomp::R
end

BooleanPartition{T}(part::R, decomp::R) where {T,R} =
    BooleanPartition{T,R}(part, decomp)

(*)(D::BooleanPartition, x::AbstractVector) = x
# reshape(x, D.decomp)

