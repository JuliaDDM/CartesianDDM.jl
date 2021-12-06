
laplacian(n::Vararg{Int}) = laplacian(Dirichlet, n...)

function laplacian(::Type{Dirichlet}, n::Int)
    spdiagm(-1 => -ones(n-1),
            0 => 2ones(n),
            1 => -ones(n-1))
end

function laplacian(::Type{Dirichlet}, n::Vararg{Int,2})
    ntot = prod(n)

    diag = view(repeat([i ≠ 1 for i in 1:n[1]], n[2]), 2:ntot)

    spdiagm(-n[1] => -ones(ntot-n[1]),
            -1 => -diag,
            0 => 4ones(ntot),
            1 => -diag,
            n[1] => -ones(ntot-n[1]))
end

