"""
    CartesianDecomposition(indices, nover, nproc)

# Examples

```julia-repl
julia> CartesianDecomposition((2:4, 3:8), (2, 1), (2, 3))
2×3 CartesianDecomposition{2, Int64, Tuple{UnitRange{Int64}, UnitRange{Int64}}}:
 (2:3, 3:5)  (2:3, 5:7)  (2:3, 7:8)
 (2:4, 3:5)  (2:4, 5:7)  (2:4, 7:8)
```

"""
struct CartesianDecomposition{N,T,R<:NTuple{N,AbstractUnitRange{T}}} <: AbstractArray{R,N}
    indices::R
    nover::NTuple{N,T}
    nproc::NTuple{N,T}
end

getindices(p::CartesianDecomposition) = p.indices
getnover(p::CartesianDecomposition) = p.nover
getnproc(p::CartesianDecomposition) = p.nproc

size(p::CartesianDecomposition) = getnproc(p)

getindex(part::CartesianDecomposition, index::CartesianIndex) =
    getindex(part, Tuple(index)...) 
function getindex(part::CartesianDecomposition, index...)
    indices = getindices(part)
    nover = getnover(part)
    nproc = getnproc(part)

    ls = first.(indices)
    hs = last.(indices)
    Ts = typeof.(indices)

    map(Ts, ls, hs, nover, nproc, index) do T, l, h, o, p, i
        d = h - l + 1
        T(
          max(d * (i - 1) ÷ p + l - o ÷ 2, l),
          min(d * i ÷ p - 1 + l + (o +      1) ÷ 2, h)
         )
    end
end


