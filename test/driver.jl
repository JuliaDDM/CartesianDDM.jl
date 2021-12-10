using CartesianDDM

dofs = (7, 12)
procs = (3, 4)
overs = (1, 1)

dims = CartesianDDMDimension.(dofs, procs, overs)

glb = globalranges(dims)
cnt = continuousranges(dims)
lcl = localranges(dims)

glb, cnt, lcl

