module Thomas

using SparseArrays

import Base: size, getindex

export laplacian, CartesianDecomposition

include("boundarycondition.jl")
include("laplacian.jl")
include("cartesian.jl")

end
