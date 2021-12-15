module CartesianDDM

using Base.Iterators
using LinearAlgebra
using OffsetArrays

import Base: OneTo
import Base: parent, size, getindex, setindex!, similar
import Base: *, \

import LinearAlgebra: ldiv!

export CartesianDDMDimension
export CartesianDDMContext
export CartesianDDMVector
export CartesianDDMRAS

export decompose

include("decomposition.jl")
include("arrays.jl")
include("vectors.jl")
include("matrices.jl")
include("coherence.jl")
include("preconditioner.jl")
include("ras.jl")

end
