module CartesianDDM

using Base.Iterators
using SparseArrays
using OffsetArrays

import Base: size, getindex

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
