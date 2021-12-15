using LinearAlgebra
using SparseArrays
using IterativeSolvers
using CartesianCutCell

include("driver.jl")

A = laplacian(dofs...)

b = rand(prod(dofs))
x = zeros(prod(dofs))

db = decompose(context, b)
dx = decompose(context, x)

cg!(x, A, b)
cg!(dx, A, db)

@assert iszero(x - dx)

dA = decompose(context, A)

@assert iszero(dA * db - A * b)

dc = CartesianDDMVector(rand, context)
c = collect(dc)

@assert iszero(dA * dc - A * c)

