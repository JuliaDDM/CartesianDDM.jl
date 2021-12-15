using CartesianDDM

dofs = (8, 12)
procs = (2, 3)
overs = (3, 2)

dims = CartesianDDMDimension.(dofs, procs, overs)
context = CartesianDDMContext(dims)

v = rand(prod(dofs))
dv = decompose(context, v)

dv[3] = 2
w = collect(dv)

