using CartesianDDM

dofs = (7, 12)
procs = (3, 4)
overs = (1, 1)

dims = CartesianDDMDimension.(dofs, procs, overs)
context = CartesianDDMContext(dims)

u = CartesianDDMVector(zeros, context)
u[12] = 1

nothing
