using CartesianDDM

dofs = (8,)
procs = (2,)
overs = (3,)

dims = CartesianDDMDimension.(dofs, procs, overs)
context = CartesianDDMContext(dims)

v = rand(prod(dofs))
dv = decompose(context, v)

dv[3] = 2
w = collect(dv)

