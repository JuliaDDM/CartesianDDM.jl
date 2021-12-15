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

@assert isapprox(x, dx)

dA = decompose(context, A)

@assert isapprox(dA * db, A * b)

dc = CartesianDDMVector(rand, context)
c = collect(dc)

@assert isapprox(dA * dc, A * c)

dP = CartesianDDMRAS(dA)

x .= 0
dx = decompose(context, x)

# bicgstabl! : does not work
# gmres!
# cg! : problem is that RAS is not symmetric
cg!(dx, A, db, Pl = dP)
cg!(x, A, b)

@show norm(x - dx)

# need to write own stationary method (block jacobi)

