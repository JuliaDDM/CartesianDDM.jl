using Thomas

indices = (1:4, 1:8)
nover = (1, 1)
nproc = (2, 3)

part = collect(CartesianDecomposition(indices,
                                      ntuple(zero, length(nover)),
                                      nproc))
decomp = collect(CartesianDecomposition(indices,
                                        nover,
                                        nproc))

