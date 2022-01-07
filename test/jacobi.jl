using CartesianDDM
using Printf
using LinearAlgebra
using SparseArrays
using IterativeSolvers
using CartesianCutCell

dofs = (8, 12)
procs = (2, 3)
overs = (4, 2)

dims = CartesianDDMDimension.(dofs, procs, overs)
context = CartesianDDMContext(dims)

A = laplacian(dofs)
dA = decompose(context, A)

P = Diagonal(A)
dP = CartesianDDMRAS(dA)

rhs = rand(prod(dofs))
drhs = decompose(context, rhs)

u = zeros(prod(dofs))
du = decompose(context, u)

res = similar(rhs)
dres = similar(drhs)

io = stdout

str = "%12s%2s" * "%12s" ^ 2 * "\n"
fmt = Printf.Format(str)

Printf.format(io, fmt, "Step", "  ", "Diagonal", "Block diag.")

str = "%12d%2s" * "%12.5E" ^ 2 * "\n"
fmt = Printf.Format(str)

itmax = 50

for it in 1:itmax
    res .= rhs - A * u
    dres .= drhs - dA * du

    u .+= P \ res
    du .+= dP \ dres

    Printf.format(io, fmt, it, "  ", norm(res), norm(dres))
end

