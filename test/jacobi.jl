using CartesianDDM
using LinearAlgebra
using SparseArrays
using IterativeSolvers
using CartesianCutCell

dofs = (8, 12)
procs = (2, 3)
overs = (4, 2)

dims = CartesianDDMDimension.(dofs, procs, overs)
context = CartesianDDMContext(dims)

A = laplacian(dofs...)
dA = decompose(context, A)

P = Diagonal(A)
dP = CartesianDDMRAS(dA)

rhs = rand(prod(dofs))
drhs = decompose(context, rhs)

u = zeros(prod(dofs))
du = decompose(context, u)

res = similar(rhs)
dres = similar(drhs)

itmax = 50

for it in 1:itmax
    res .= rhs - A * u
    dres .= drhs - dA * du

    u .+= P \ res
    du .+= dP \ dres

    println("it=$it : $(norm(res)) vs. $(norm(dres))")
end

