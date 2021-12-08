module CartesianDDM

using Base.Iterators
using OffsetArrays

import Base: size, getindex, setindex!

export decompose, partition
export CartesianDDMVector
export CartesianBooleanPartition

include("decomposition.jl")
include("arrays.jl")
include("vectors.jl")
include("matrices.jl")
include("partition.jl")
include("synchronize.jl")

end
