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

