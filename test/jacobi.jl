### A Pluto.jl notebook ###
# v0.17.2

using Markdown
using InteractiveUtils

# ╔═╡ cb9498e0-5835-11ec-0956-51f67d2e3f79
begin
    import Pkg

    # activate the shared project environment
    Pkg.activate(Base.current_project())
    # instantiate, i.e. make sure that all packages are downloaded
    Pkg.instantiate()
    using LinearAlgebra, SparseArrays, CartesianDDM

	Pkg.add("IterativeSolvers")
	using IterativeSolvers

	Pkg.add(url="https://github.com/JuliaCutCell/CartesianCutCell.jl.git")
	using CartesianCutCell
end


# ╔═╡ db9b6205-f8cc-443c-b5d9-a4d6ee74ff05
begin
	indices = (3:4, 2:3)
	nproc = (3, 4)
	nover = (1, 1)
	indices, nproc, nover
end

# ╔═╡ bf3a023b-3966-4ed5-812e-7249dc8e36d9
begin
	decomp = decompose(indices, nproc, nover)
	part = decompose(indices, nproc)
	decomp, part
end

# ╔═╡ 55d3914f-f780-4ad1-92f5-1647d7534b01
begin
	b = CartesianDDMVector(rand, decomp, part)
	x = CartesianDDMVector(zeros, decomp, part)
end

# ╔═╡ 53d98cef-4e56-42ae-8b7e-b67dbd1632d4
A = laplacian(length.(indices))

# ╔═╡ 3dc10aee-e5cf-486f-bc7c-67f4802e5f8d
Pl = Diagonal(A)

# ╔═╡ 8247be53-13b2-4dd1-b281-dcdb275a5bb3
cg!(x, A, b; Pl)

# ╔═╡ Cell order:
# ╠═cb9498e0-5835-11ec-0956-51f67d2e3f79
# ╠═db9b6205-f8cc-443c-b5d9-a4d6ee74ff05
# ╠═bf3a023b-3966-4ed5-812e-7249dc8e36d9
# ╠═55d3914f-f780-4ad1-92f5-1647d7534b01
# ╠═53d98cef-4e56-42ae-8b7e-b67dbd1632d4
# ╠═3dc10aee-e5cf-486f-bc7c-67f4802e5f8d
# ╠═8247be53-13b2-4dd1-b281-dcdb275a5bb3
