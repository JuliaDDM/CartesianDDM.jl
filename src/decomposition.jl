"""
    abstract type AbstractIndexing end

Abstract indexing context.

"""
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

- `dof` : number of degrees of freedom,
- `proc` : number of domains,
- `over` : overlap.

# To do list

- Include periodicity.

"""
struct CartesianDDMDimension
    dof::Int
    proc::Int
    over::Int
end

getdof(dim::CartesianDDMDimension) = getproperty(dim, :dof)
getproc(dim::CartesianDDMDimension) = getproperty(dim, :proc)
getover(dim::CartesianDDMDimension) = getproperty(dim, :over)

function removeoverlap(dim::CartesianDDMDimension)
    (; dof, proc, over) = dim
    CartesianDDMDimension(dof, proc, 0)
end

struct CartesianDDMContext{N}# <: StaticVector{N,CartesianDDMDimension}
    dims::NTuple{N,CartesianDDMDimension}
end

getdims(context::CartesianDDMContext) = getproperty(context, :dims)

function removeoverlap(context::CartesianDDMContext)
    (; dims) = context
    CartesianDDMContext(removeoverlap.(dims))
end

"""
    CartesianDDMRanges(indexing, dims)

Multidimensional array storing each domains' ranges according
to an indexing context.

"""
struct CartesianDDMRanges{R<:AbstractUnitRange{Int},N,I} <: AbstractArray{NTuple{N,R},N}
    indexing::I
    context::CartesianDDMContext{N}
end

getindexing(iter::CartesianDDMRanges) = getproperty(iter, :indexing)
getcontext(iter::CartesianDDMRanges) = getproperty(iter, :context)
getdims(iter::CartesianDDMRanges) = getproperty(getcontext(iter), :dims)

cartesianddmranges(indexing, context::CartesianDDMContext{N}) where {N} =
    CartesianDDMRanges{_rangetype(indexing),N,typeof(indexing)}(indexing, context)

globalranges(context) =
    cartesianddmranges(globally, context)
continuousranges(context) =
    cartesianddmranges(continuously, context)
localranges(context) =
    cartesianddmranges(locally, context)

ranges(context) =
    globalranges(context),
    continuousranges(context),
    localranges(context)

function size(iter::CartesianDDMRanges)
    context = getcontext(iter)
    dims = getdims(context)
    getproc.(dims)
end

getindex(iter::CartesianDDMRanges, index::CartesianIndex) =
    getindex(iter, Tuple(index)...)
getindex(iter::CartesianDDMRanges, index::Int...) =
    _getindex(getindexing(iter), getcontext(iter), index...)

function _getindex(::Globally, context, index...)
    dims = getdims(context)

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

function _getindex(::Continuously, context, index...) where {R}
    glbs = _getindex(globally, context, index...)

    dims = getdims(context)

    overs = getover.(dims)

    map(overs, index, glbs) do o, i, g
        _rangetype(continuously)(
          o * (i - 1) + first(g),
          o * (i - 1) + last(g)
         )
    end
end

function _getindex(::Locally, context, index...) where {R}
    glbs = _getindex(globally, context, index...)
    _rangetype(locally).(length.(glbs))
end

"""

Full range (works only for global and continuous).

"""
function span end

