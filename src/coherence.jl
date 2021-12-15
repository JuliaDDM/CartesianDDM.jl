ArrayAbstract{N,T} = AbstractArray{T,N}

"""
Assumed PoU is Boolean.

"""
makecoherent

"""
    makecoherent!(x::CartesianDDMVector)

Make coherent.

"""
function makecoherent!(x::CartesianDDMVector)
    (; context, parent) = x
    nocontext = removeoverlap(context)

    dec, _ = ranges(context)
    part, _ = ranges(nocontext)

    makecoherent!(dec, part, parent)

    x
end

"""
Make coherent in 1d

"""
#=
makecoherent!(dec::ArrayAbstract{1},
              part::ArrayAbstract{1},
              parent::ArrayAbstract{1}) =
    _makecoherent!(dec, part, parent)
=#

"""
Make coherent across faces

"""
function makecoherent!(dec::ArrayAbstract{1},
                       part::ArrayAbstract{1},
                       parent::ArrayAbstract{1})
    for i in drop(axes(parent, 1), 1)
        _makecoherent!((parent[i-1], dec[i-1]),
                       (parent[i], dec[i], part[i]))
        _makecoherent!((parent[i], dec[i]),
                       (parent[i-1], dec[i-1], part[i-1]))
    end
end

"""
Make coherent in 2d

"""
function makecoherent!(dec::ArrayAbstract{2},
                       part::ArrayAbstract{2},
                       parent::ArrayAbstract{2})
    for j in axes(parent, 2)
        makecoherent!(@view(dec[:, j]),
                      @view(part[:, j]),
                      @view(parent[:, j]))
    end
    for i in axes(parent, 1)
        makecoherent!(@view(dec[i, :]),
                      @view(part[i, :]),
                      @view(parent[i, :]))
    end

    _makecoherent!(dec, part, parent)
end

"""
Make coherent across edges

"""
function _makecoherent!(dec::ArrayAbstract{2},
                        part::ArrayAbstract{2},
                        parent::ArrayAbstract{2})
    for j in drop(axes(parent, 2), 1)
        for i in drop(axes(parent, 1), 1)
            _makecoherent!((parent[i-1, j-1], dec[i-1, j-1]),
                           (parent[i, j], dec[i, j], part[i, j]))
            _makecoherent!((parent[i, j], dec[i, j]),
                           (parent[i-1, j-1], dec[i-1, j-1], part[i-1, j-1]))
            _makecoherent!((parent[i, j-1], dec[i, j-1]),
                           (parent[i-1, j], dec[i-1, j], part[i-1, j]))
            _makecoherent!((parent[i-1, j], dec[i-1, j]),
                           (parent[i, j-1], dec[i, j-1], part[i, j-1]))
        end
    end
end

"""
Make coherent in 3d

"""
function makecoherent!(dec::ArrayAbstract{3},
                       part::ArrayAbstract{3},
                       parent::ArrayAbstract{3})
    for k in axes(parent, 3)
        for j in axes(parent, 2)
            makecoherent!(@view(dec[:, j, k]),
                          @view(part[:, j, k]),
                          @view(parent[:, j, k]))
        end
    end
    for k in axes(parent, 3)
        for i in axes(parent, 1)
            makecoherent!(@view(dec[i, :, k]),
                          @view(part[i, :, k]),
                          @view(parent[i, :, k]))
        end
    end
    for j in axes(parent, 2)
        for i in axes(parent, 1)
            makecoherent!(@view(dec[i, j, :]),
                          @view(part[i, j, :]),
                          @view(parent[i, j, :]))
        end
    end

    for k in axes(parent, 3)
        makecoherent!(@view(dec[:, :, k]),
                      @view(part[:, :, k]),
                      @view(parent[:, :, k]))
    end
    for j in axes(parent, 2)
        makecoherent!(@view(dec[:, j, :]),
                      @view(part[:, j, :]),
                      @view(parent[:, j, :]))
    end
    for i in axes(parent, 1)
        makecoherent!(@view(dec[i, :, :]),
                      @view(part[i, :, :]),
                      @view(parent[i, :, :]))
    end

    _makecoherent!(dec, part, parent)
end

"""
Make coherent across vertices

"""
function _makecoherent!(dec::ArrayAbstract{3},
                        part::ArrayAbstract{3},
                        parent::ArrayAbstract{3})
    for k in drop(axes(parent, 3), 1)
        for j in drop(axes(parent, 2), 1)
            for i in drop(axes(parent, 1), 1)
                _makecoherent!((parent[i-1, j-1, k-1], dec[i-1, j-1, k-1]),
                               (parent[i, j, k], dec[i, j, k], part[i, j, k]))
                _makecoherent!((parent[i, j, k], dec[i, j, k]),
                               (parent[i-1, j-1, k-1], dec[i-1, j-1, k-1], part[i-1, j-1, k-1]))
                _makecoherent!((parent[i, j, k-1], dec[i, j, k-1]),
                               (parent[i-1, j-1, k], dec[i-1, j-1, k], part[i-1, j-1, k]))
                _makecoherent!((parent[i-1, j-1, k], dec[i-1, j-1, k]),
                               (parent[i, j, k-1], dec[i, j, k-1], part[i, j, k-1]))
                _makecoherent!((parent[i-1, j, k-1], dec[i-1, j, k-1]),
                               (parent[i, j-1, k], dec[i, j-1, k], part[i, j-1, k]))
                _makecoherent!((parent[i, j-1, k], dec[i, j-1, k]),
                               (parent[i-1, j, k-1], dec[i-1, j, k-1], part[i-1, j, k-1]))
                _makecoherent!((parent[i, j-1, k-1], dec[i, j-1, k-1]),
                               (parent[i-1, j, k], dec[i-1, j, k], part[i-1, j, k]))
                _makecoherent!((parent[i-1, j, k], dec[i-1, j, k]),
                               (parent[i, j-1, k-1], dec[i, j-1, k-1], part[i, j-1, k-1]))
            end
        end
    end
end

"""

"""
function _makecoherent!((out, gblout), (in, gblin, nogblin))
    rout = reshape(out, gblout...)
    rin = reshape(in, gblin...)

    iter = intersect(CartesianIndices.((gblout, nogblin))...)

    rout[iter] .= rin[iter]
end

