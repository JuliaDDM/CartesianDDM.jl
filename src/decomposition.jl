abstract type AbstractIndexing end

struct Globally <: AbstractIndexing end
struct Continuously <: AbstractIndexing end
struct Locally <: AbstractIndexing end

const globally = Globally()
const continuously = Continuously()
const locally = Locally()

_rangetype(::Globally) = UnitRange{Int}
_rangetype(::Continuously) = UnitRange{Int}
_rangetype(::Locally) = Base.OneTo{Int}

"""

    CartesianDDMDimension(dof, proc, over)

Stores directional information.

# To do list

1. Include periodicity.

"""
struct CartesianDDMDimension
    dof::Int
    proc::Int
    over::Int
end

getdof(dim::CartesianDDMDimension) = getproperty(dim, :dof)
getproc(dim::CartesianDDMDimension) = getproperty(dim, :proc)
getover(dim::CartesianDDMDimension) = getproperty(dim, :over)

#=
"""
    CartesianDecomposition(ranges, nproc, nover)

Aims at storing global indices of each subdomain.

Is it useful to allow `ranges` to be something else than `NTuple{N,Base.OneTo}`?
If not, just replace `ranges` by `ndof` and be done with it.

"""
struct CartesianDecomposition{N,I,R<:NTuple{N,AbstractUnitRange{Int}}} <: AbstractArray{R,N}
    indexing::I
    ranges::R
    nproc::NTuple{N,Int}
    nover::NTuple{N,Int}
end
=#

struct CartesianDDMRanges{R<:AbstractUnitRange{Int},N,I} <: AbstractArray{NTuple{N,R},N}
    indexing::I
    dims::NTuple{N,CartesianDDMDimension}
end

getindexing(iter::CartesianDDMRanges) = getproperty(iter, :indexing)
getdims(iter::CartesianDDMRanges) = getproperty(iter, :dims)

cartesianddmranges(indexing, dims::NTuple{N}) where {N} =
    CartesianDDMRanges{_rangetype(indexing),N,typeof(indexing)}(indexing, dims)
globalranges(dims) =
    cartesianddmranges(globally, dims)
continuousranges(dims) =
    cartesianddmranges(continuously, dims)
localranges(dims) =
    cartesianddmranges(locally, dims)

size(iter::CartesianDDMRanges) = getproc.(getdims(iter))

getindex(iter::CartesianDDMRanges, index::CartesianIndex) =
    getindex(iter, Tuple(index)...)

getindex(iter::CartesianDDMRanges, index::Int...) =
    _getindex(getindexing(iter), getdims(iter), index...)

function _getindex(::Globally, dims, index...) where {R}
    dofs = getdof.(dims)
    procs = getproc.(dims)
    overs = getover.(dims)

    map(dofs, procs, overs, index) do d, p, o, i
        _rangetype(globally)(
          max(d * (i - 1) ÷ p + 1 - o ÷ 2, 1),
          min(d * i ÷ p - 1 + 1 + (o + 1) ÷ 2, d)
         )
    end
end

function _getindex(::Continuously, dims, index...) where {R}
    glbs = _getindex(globally, dims, index...)

    overs = getover.(dims)

    map(overs, index, glbs) do o, i, g
        _rangetype(continuously)(
          o * (i - 1) + first(g),
          o * (i - 1) + last(g)
         )
    end
end

function _getindex(::Locally, dims, index...) where {R}
    glbs = _getindex(globally, dims, index...)
    _rangetype(locally).(length.(glbs))
end

"""

Full range (works only for global and continuous).

"""
function span end

#=
decompose(ranges, nproc) =
    decompose(ranges, nproc, ntuple(zero, length(ranges)))

decompose(ranges, nproc, nover, indexing=globally) =
    CartesianDecomposition(indexing, ranges, nproc, nover)

size(iter::CartesianDecomposition) = getproperty(iter, :nproc)

getindex(iter::CartesianDecomposition, index::CartesianIndex) =
    getindex(iter, Tuple(index)...)

function getindex(iter::CartesianDecomposition{N}, index::Vararg{Int,N}) where {N}
    (; indexing, ranges, nproc, nover) = iter
    _getindex(indexing, ranges, nproc, nover, index...)
end

function _getindex(::Globally, ranges, nproc, nover, index...)
    ls = first.(ranges)
    hs = last.(ranges)
    Ts = typeof.(ranges)

    map(Ts, ls, hs, nover, nproc, index) do T, l, h, o, p, i
        d = h - l + 1
        T(
          max(d * (i - 1) ÷ p + l - o ÷ 2, l),
          min(d * i ÷ p - 1 + l + (o + 1) ÷ 2, h)
         )
    end
end

"""
Should it start from ``1``?

"""
function _getindex(::Continuously, ranges, nproc, nover, index...)
    el = _getindex(globally, ranges, nproc, nover, index...)
    map(nover, index, el) do o, i, r
        typeof(r)(o * (i - 1) + first(r),
                  o * (i - 1) + last(r))
    end
end

"""
Do not use `axes.(el, 1)` to guarantee `eltype(iter) == R`.
Work on `CartesianDecomposition` to fix this.

"""
function _getindex(::Locally, ranges, nproc, nover, index...)
    el = _getindex(globally, ranges, nproc, nover, index...)
    map(el) do r
        typeof(r)(1, length(r))
    end
end

indexglobally(iter) = _switchto(globally, iter)
indexcontinuously(iter) = _switchto(continuously, iter)
indexlocally(iter) = _switchto(locally, iter)

function _switchto(this, iter)
    (; indexing, ranges, nproc, nover) = iter
    CartesianDecomposition(this, ranges, nproc, nover)
end
=#

