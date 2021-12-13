module CartesianDDM

using Base.Iterators
using OffsetArrays

import Base: OneTo, parent, size, getindex, setindex!

export CartesianDDMDimension
export CartesianDDMContext
export CartesianDDMVector

include("decomposition.jl")
include("arrays.jl")
include("vectors.jl")

#=
include("matrices.jl")
include("partition.jl")
include("synchronize.jl")
=#

end
