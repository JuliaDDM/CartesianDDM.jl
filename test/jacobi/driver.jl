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

# do something inside of U̲
synchronize!((loc, dec), (glob, part))

nline = 2(3prod(nproc)
          - nproc[1] * nproc[2]
          - nproc[2] * nproc[3]
          - nproc[3] * nproc[1])
nedge = 4(3prod(nproc)
          -2(nproc[1] * nproc[2]
             + nproc[2] * nproc[3]
             + nproc[3] * nproc[1])
          + sum(nproc))
nvertex = 8prod(nproc .- 1)

@show nline, nedge, nvertex

