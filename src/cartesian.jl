"""
    CartesianStarMatrix{S,N,T,AA} <: AbstractMatrix{T}

"""
struct CartesianStarMatrix{S,N,T,AA<:NTuple{S,NTuple{N,AbstractArray{T,N}}}} <: AbstractMatrix{T}
    diags::AA
end

Base.size(A::CartesianStarMatrix) = ntuple(i -> prod(size(first(first(diagonals(A))))), 2)

function Base.getindex(A::CartesianStarMatrix, i, j)
    indices = CartesianIndices(first(first(diagonals(A))))

    I, J = indices[i].I, indices[j].I

    sum(abs.(I .- J)) > 1 ? zero(eltype(A)) : one(eltype(A))
end

diagonals(A::CartesianStarMatrix) = A.diags

"""
    dirichletlaplacian(n...)

Second order laplacian with Dirichlet boundary conditions.

"""
function dirichletlaplacian(n...)

    coefs = (-1.0, 2.0, -1.0)

    diags = map(coefs) do c
        ntuple(length(n)) do i
            c * ones(n)
        end
    end

    CartesianStarMatrix(diags)
end

"""
    tridiagonal(A)

Convert a `CartesianStarMatrix{3,1}` matrix into a `LinearAlgebra.Tridiagonal` matrix.

"""
function tridiagonal(A::CartesianStarMatrix{3,1})
    ((α,), (β,), (γ,)) = diagonals(A)

    a = @view α[begin + 1:end]
    b = @view β[begin:end]
    c = @view γ[begin:end - 1]

    Tridiagonal(a, b, c)
end

*(A::CartesianStarMatrix{3,1}, x::AbstractVector) = tridiagonal(A) * x

#=
A = dirichletlaplacian(3, 4);
b = rand(12);
A * b
=#

