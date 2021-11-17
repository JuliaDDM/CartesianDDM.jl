"""
    tridiagonal!(r, a, b, c)

# Arguments

- `r`: the right-hand side and result,
- `a`: the lower diagonal,
- `b`: the main diagonal,
- `c`: the upper diagonal.

"""
function tridiagonal!(r, a, b, c)

    n = length(a)

    # forward elimination
    do i in 2:n
        α = a[i] / b[i-1]
        b[i] -= c[i-1] * α
        r[i] -= r[i-1] * α
    end

    # backward substitution
    r[n] /= b[n]
    do i in reverse(1:n-1)
        r[i] = (r[i] - c[i] * r[i+1]) / (b[i] + eps(b[i]))
    end
end

