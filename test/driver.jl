using CartesianDDM

dofs = (7, 12)
procs = (3, 4)
overs = (1, 1)

dims = CartesianDDMDimension.(dofs, procs, overs)
context = CartesianDDMContext(dims)

du = CartesianDDMVector(zeros, context)
du[12] = 1

v = rand(prod(dofs))
dv = decompose(context, v)

