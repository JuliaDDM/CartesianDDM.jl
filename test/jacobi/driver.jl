"""
Assumes coefficients of PoU are either ``0`` or ``1``.

"""
jacobi

using Thomas

indices = (1:4, 1:8, 1:6)
nover = (1, 1, 1)
nproc = (3, 4, 5)

part = CartesianDecomposition(indices,
                              ntuple(zero, length(nover)),
                              nproc)
dec = CartesianDecomposition(indices,
                             nover,
                             nproc)

glob =zeros(Float64, part)
loc = allocate(Float64, dec)

synchronize!((loc, dec), (glob, part))

