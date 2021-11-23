"""
Assumes coefficients of PoU are either ``0`` or ``1``.

"""
jacobi

using Thomas

indices = (1:4, 1:8)
nover = (1, 1)
nproc = (2, 3)

part = CartesianDecomposition(indices,
                              ntuple(zero, length(nover)),
                              nproc)
dec = CartesianDecomposition(indices,
                             nover,
                             nproc)

U̲ =zeros(Float64, dec) 
V̲ = allocate(Float64, dec)

# do something inside of U̲
synchronize!(V̲, U̲, part, dec)

