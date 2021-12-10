module CartesianDDM

using Base.Iterators
using OffsetArrays

import Base: OneTo, parent, size, getindex#, setindex!

export WrappedDDMVector
export indexglobally, indexcontinuously, indexlocally

export decompose
export CartesianDecomposition
export CartesianDDMVector
export CartesianBooleanPartition

export CartesianDDMDimension, CartesianDDMRanges
export Globally, Continuously, Locally
export globally, continuously, locally
export globalranges, continuousranges, localranges

include("decomposition.jl")
include("arrays.jl")
include("vectors.jl")
include("matrices.jl")
include("partition.jl")
include("synchronize.jl")

end
