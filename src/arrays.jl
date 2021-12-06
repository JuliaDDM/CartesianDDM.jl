function zeros(::Type{T}, iter) where {T}
    map(iter) do el
        zeros(prod(length.(el)))
    end
end

function allocate(::Type{T}, iter) where {T}
    map(iter) do el
        Vector{T}(undef, prod(length.(el)))
    end
end

#=
struct CartesianVector{T,R,A<:AbstractVector{T}} <: AbstractVector{T}
    indices::R
    parent::A
end

function CartesianVector{T}(::UndefInitializer, indices) where {T}
    parent = Vector{T}(undef, prod(length.(indices)))
    CartesianVector(indices, parent)
end

function CartesianVector(init::Function, indices)
    parent = init(prod(length.(indices)))
    CartesianVector(indices, parent)
end

parent(x::CartesianVector) = getproperty(x, :parent)
size(x::CartesianVector) = size(parent(x))
getindex(x::CartesianVector, i...) = getindex(parent(x), i...)
(*)(A::AbstractMatrix, x::CartesianVector) = A * parent(x)

=#

