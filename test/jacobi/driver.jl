"""
Assumes coefficients of PoU are either ``0`` or ``1``.

"""
jacobi

using CartesianDDM

indices = (1:4, 1:8, 1:6)
nover = (1, 1, 1)
nproc = (3, 4, 5)

decomp = decompose(indices, nover, nproc)
Ui = CartesianDDMVector(zeros, decomp)

part = partition(indices, nproc)
Di = CartesianBooleanPartition{Float64}(part)

nothing
#Di .* Ui

#=
glob =zeros(Float64, part)
loc = allocate(Float64, decomp)

synchronize!((loc, decomp), (glob, part))
=#

