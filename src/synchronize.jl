ArrayAbstract{N,T} = AbstractArray{T,N}

"""
Synchronize in 1d

"""
function synchronize!((out, dec)::T,
                      (in, part)::T) where {T<:NTuple{2,ArrayAbstract{1}}}
    _synchronize!((out, dec), (in, part))
end

"""
Synchronize in 2d

"""
function synchronize!((out, dec)::T,
                      (in, part)::T) where {T<:NTuple{2,ArrayAbstract{2}}}
    for j in axes(out, 2)
        _synchronize!((@view(out[:, j]), @view(dec[:, j])),
                      (@view(in[:, j]), @view(part[:, j])))
    end
    for i in axes(out, 1)
        _synchronize!((@view(out[i, :]), @view(dec[i, :])),
                      (@view(in[i, :]), @view(part[i, :])))
    end

    _synchronize!((out, dec), (in, part))
end

"""
Synchronize in 3d

"""
function synchronize!((out, dec)::T,
                      (in, part)::T) where {T<:NTuple{2,ArrayAbstract{3}}}
    for k in axes(out, 3)
        for j in axes(out, 2)
            _synchronize!((@view(out[:, j, k]), @view(dec[:, j, k])),
                          (@view(in[:, j, k]), @view(part[:, j, k])))
        end
    end
    for k in axes(out, 3)
        for i in axes(out, 1)
            _synchronize!((@view(out[i, :, k]), @view(dec[i, :, k])),
                          (@view(in[i, :, k]), @view(part[i, :, k])))
        end
    end
    for j in axes(out, 2)
        for i in axes(out, 1)
            _synchronize!((@view(out[i, j, :]), @view(dec[i, j, :])),
                          (@view(in[i, j, :]), @view(part[i, j, :])))
        end
    end

    for k in axes(out, 3)
        _synchronize!((@view(out[:, :, k]), @view(dec[:, :, k])),
                      (@view(in[:, :, k]), @view(part[:, :, k])))
    end
    for j in axes(out, 2)
        _synchronize!((@view(out[:, j, :]), @view(dec[:, j, :])),
                      (@view(in[:, j, :]), @view(part[:, j, :])))
    end
    for i in axes(out, 1)
        _synchronize!((@view(out[i, :, :]), @view(dec[i, :, :])),
                      (@view(in[i, :, :]), @view(part[i, :, :])))
    end

    _synchronize!((out, dec), (in, part))
end

"""
Synchronize across faces

"""
function _synchronize!((out, dec)::T,
                       (in, part)::T) where {T<:NTuple{2,ArrayAbstract{1}}}
    for i in drop(axes(out, 1), 1)
        _synchronize!((out[i-1], dec[i-1]),
                      (in[i], part[i]), 1)
        _synchronize!((out[i], dec[i]),
                      (in[i-1], part[i-1]), 1)
    end
end

"""
Synchronize across edges

"""
function _synchronize!((out, dec)::T,
                       (in, part)::T) where {T<:NTuple{2,ArrayAbstract{2}}}
    for j in drop(axes(out, 2), 1)
        for i in drop(axes(out, 1), 1)
            _synchronize!((out[i-1, j-1], dec[i-1, j-1]),
                          (in[i, j], part[i, j]), 2)
            _synchronize!((out[i, j], dec[i, j]),
                          (in[i-1, j-1], part[i-1, j-1]), 2)
            _synchronize!((out[i, j-1], dec[i, j-1]),
                          (in[i-1, j], part[i-1, j]), 2)
            _synchronize!((out[i-1, j], dec[i-1, j]),
                          (in[i, j-1], part[i, j-1]), 2)
        end
    end
end

"""
Synchronize across vertices

"""
function _synchronize!((out, dec)::T,
                       (in, part)::T) where {T<:NTuple{2,ArrayAbstract{3}}}
    for k in drop(axes(out, 3), 1)
        for j in drop(axes(out, 2), 1)
            for i in drop(axes(out, 1), 1)
                _synchronize!((out[i-1, j-1, k-1], dec[i-1, j-1, k-1]),
                              (in[i, j, k], part[i, j, k]), 3)
                _synchronize!((out[i, j, k], dec[i, j, k]),
                              (in[i-1, j-1, k-1], part[i-1, j-1, k-1]), 3)
                _synchronize!((out[i, j, k-1], dec[i, j, k-1]),
                              (in[i-1, j-1, k], part[i-1, j-1, k]), 3)
                _synchronize!((out[i-1, j-1, k], dec[i-1, j-1, k]),
                              (in[i, j, k-1], part[i, j, k-1]), 3)
                _synchronize!((out[i-1, j, k-1], dec[i-1, j, k-1]),
                              (in[i, j-1, k], part[i, j-1, k]), 3)
                _synchronize!((out[i, j-1, k], dec[i, j-1, k]),
                              (in[i-1, j, k-1], part[i-1, j, k-1]), 3)
                _synchronize!((out[i, j-1, k-1], dec[i, j-1, k-1]),
                              (in[i-1, j, k], part[i-1, j, k]), 3)
                _synchronize!((out[i-1, j, k], dec[i-1, j, k]),
                              (in[i, j-1, k-1], part[i, j-1, k-1]), 3)
            end
        end
    end
end

"""
Synchronize a local array from a global array

"""
function _synchronize!((out, dec), (in, part), arg)
    cnt[arg] += 1
    println(cnt)
end

