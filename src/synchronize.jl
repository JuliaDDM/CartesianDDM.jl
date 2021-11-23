@generated function synchronize!(buf, arg...)
    N = ndims(buf)
    quote
        @nexprs $N d -> synchronize!(Val{d}, buf, arg...)
    end
end

@generated function synchronize!(::Type{Val{D}}, out, in, part, dec) where {D}
    N = ndims(out)
    quote
        @nloops $N i d -> drop(axes(out, d), d == $D) begin
            _synchronize!(@nref($N, out, d -> d == $D ? i_d-1 : i_d),
                          @nref($N, in, i),
                          @nref($N, part, i),
                          @nref($N, dec, d -> d == $D ? i_d-1 : i_d))
            _synchronize!(@nref($N, out, i),
                          @nref($N, in, d -> d == $D ? i_d-1 : i_d),
                          @nref($N, part, d -> d == $D ? i_d-1 : i_d),
                          @nref($N, dec, i))
        end
    end
end

function _synchronize!(out, in, part, dec)
    println("coucou")
end

