module CartesianDDM

using Base.Cartesian
using Base.Iterators
using SparseArrays
using OffsetArrays

import Base: zeros# delete
import Base: size, getindex, *
#import Base: parent

#export CartesianVector
export BooleanPartition
export laplacian, allocate, decompose, partition, synchronize!

const cnt = Vector{Int}(undef, 3)

include("boundarycondition.jl")
include("operators.jl")
include("cartesian.jl")
include("arrays.jl")
include("partition.jl")
include("synchronize.jl")

__init__() = cnt .= 0

end
