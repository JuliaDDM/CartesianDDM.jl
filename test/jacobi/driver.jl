"""
Assumes coefficients of PoU are either ``0`` or ``1``.

"""
jacobi

using CartesianDDM

indices = (3:4, 2:3)
nover = (1, 1)
nproc = (3, 4)

decomp = decompose(indices, nover, nproc)
part = partition(indices, nproc)

b = CartesianDDMVector(rand, decomp, part)
x = CartesianDDMVector(zeros, decomp, part)

#=
using CartesianCutCell

A = laplacian(length.(indices))
=#
#=
using LinearAlgebra, IterativeSolvers

Pl = Diagonal(A)
cg!(x, A, b; Pl)
=#

