const DDMVector{T} = DDMArray{T,1}

"""

Not a dense array!

"""
struct CartesianDDMVector{T,A} <: DDMVector{T}
    parent::A
end

function CartesianDDMVector{T}(::UndefInitializer, decomp) where {T}
    parent = map(decomp) do indices
        Vector{T}(undef, prod(length.(indices)))
    end
    CartesianDDMVector{T,typeof(parent)}(parent)
end

function CartesianDDMVector(init::Function, decomp)
    parent = map(decomp) do indices
        init(prod(length.(indices)))
    end
    CartesianDDMVector{eltype(first(parent)),typeof(parent)}(parent)
end

