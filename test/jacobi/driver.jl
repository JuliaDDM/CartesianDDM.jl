"""
Assumes coefficients of PoU are either ``0`` or ``1``.

"""
jacobi

using CartesianDDM

indices = (3:6, 2:9, 3:8)
nover = (1, 1, 1)
nproc = (3, 4, 5)

decomp = decompose(indices, nover, nproc)
part = partition(indices, nproc)

Ui = CartesianDDMVector{Float64}(undef, decomp, part)

#=
foo(ci::CartesianIndex) = prod(ci.I)

for (ar, indices) in zip(Ui.parent, decomp)
    cindices = CartesianIndices(indices)
    ar .= reshape(foo.(cindices), :)
end

Vi = reshape(foo.(CartesianIndices(indices)), indices)

for (ar, indices) in zip(Ui.parent, decomp)
    a = reshape(ar, length.(indices))
    b = Vi[CartesianIndices(indices)]
end

for (u, v) in zip(Ui, Vi)
    !iszero(u - v) && @show u - v
end
=#

#=
Di = CartesianBooleanPartition{Float64}(part)
=#

#=
glob =zeros(Float64, part)
loc = allocate(Float64, decomp)

synchronize!((loc, decomp), (glob, part))
=#

