struct CartesianDDMRAS{T,C,A} <: DDMPreconditioner{T}
    context::C
    parent::A
end

parent(P::CartesianDDMRAS) = getproperty(P, :parent)
getcontext(P::CartesianDDMRAS) = getproperty(P, :context)

function CartesianDDMRAS(A::CartesianDDMMatrix, init=factorize)
    (; context, parent) = A

    out = map(parent) do el
        init(el)
    end

    CartesianDDMRAS{eltype(A),typeof(context),typeof(out)}(context, out)
end

function size(x::CartesianDDMRAS)
    (; context, parent) = x
    (; dims) = context
    ntuple(i -> mapreduce(getdof, *, dims), 2)
end

"""
    ldiv!(y, P, x)

Computes `P \\ x` in-place of `y`.

"""
function ldiv!(y::CartesianDDMVector, P::CartesianDDMRAS, x::CartesianDDMVector)
    map(parent(y), parent(P), parent(x)) do y_, P_, x_
        y_ .= P_ \ x_
    end
    makecoherent!(y)
    y
end

"""
    ldiv!(P, x)

Computes `P \\ x` in-place of `x`.

"""
function ldiv!(P::CartesianDDMRAS, x::CartesianDDMVector)
    map(parent(P), parent(x)) do P_, x_
        x_ .= P_ \ x_
    end
    makecoherent!(x)
    x
end

"""
    \\(P, x)

Computes `P \\ x`.

"""
function (\)(P::CartesianDDMRAS, x::CartesianDDMVector)
    y = similar(x)
    ldiv!(y, P, x)
end

