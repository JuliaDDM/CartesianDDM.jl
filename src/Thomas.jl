module Thomas

using SparseArrays
using Base.Cartesian
using Base.Iterators

import Base: size, getindex, zeros

export laplacian, CartesianDecomposition, allocate, synchronize!

include("boundarycondition.jl")
include("laplacian.jl")
include("cartesian.jl")
include("arrays.jl")
include("synchronize.jl")

end
