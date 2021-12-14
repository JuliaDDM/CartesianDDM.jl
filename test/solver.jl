using LinearAlgebra
#using SparseArrays
using IterativeSolvers
using CartesianCutCell

include("driver.jl")

A = laplacian(dofs...)

b = rand(prod(dofs))
x = zeros(prod(dofs))

# synchronization/make coherent is missing
db = decompose(context, b)
dx = decompose(context, x)

cg!(x, A, b)
cg!(dx, A, db)

@assert iszero(x - dx)

dA = decompose(context, A)

@assert iszero(dA * db - A * b)

#=
rows = rowvals(A)
vals = nonzeros(A)
m, n = size(A)
for j = 1:n
    for i in nzrange(A, j)
        row = rows[i]
        val = vals[i]
        dval = dA[row, j]
        !iszero(val - dval) && @show val, dval
    end
end
=#

