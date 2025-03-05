### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 980691f0-3ec3-11ef-0ef1-e7653399f7de
begin
	using Pkg; Pkg.activate("../.")
	using FileIO, JLD2, TemporalNetworks
end

# ╔═╡ 7324ec91-a825-47c4-b2aa-2d8a22b5ba3c
# Example 1: Transition from unclustered to a single appearing cluster

# ╔═╡ b8449181-55d3-42e1-83cf-57ede95e1033
begin
	W = FileIO.load("Ex1.jld2","W") |> Vector{Matrix{Float64}}
	mlgraph = MultilayerGraph(W, connect = Multiplex())
	partition = SpectralPartition(mlgraph)
	seba_part = SEBAPartition(partition,1)
end

# ╔═╡ 3f621722-cced-45f8-81ce-485234c289a3
# Plot the network at its first and final time steps

# ╔═╡ 8f50edf1-fe0c-4e43-8fa6-d02639c22416
begin
	p1 = [plot(mlgraph,x) for x in [1,21]]
	plot(p1...)
end

# ╔═╡ dd9c3a1e-6bad-4895-8342-33b6fa9a43aa
# Plot leading evecs of the supra Laplacian

# ╔═╡ cb73e5c2-b133-430e-9ca3-d052b7262bba
plot(plot(partition)[1][1:6]...)

# ╔═╡ 15b20c67-a92b-401c-81ff-11e7209afdd2
# Perform SEBA (with doubling) of the leading spatial eigenvector.

# ╔═╡ daa0159f-da2b-4bb3-8abc-51ed89124b25
plot(plot(seba_part)[1]...)

# ╔═╡ 9d67cfd0-d114-4bc9-b765-afe8b2221469
# Example 2: Transition from unclustered to single appearing cluster that disintegrates into two clusters further.

# ╔═╡ 5693c8ae-c7cb-425f-91b0-5ac8d38a1c55
begin
	W_2 = FileIO.load("Ex2.jld2","W") |> Vector{Matrix{Float64}}
	mlgraph_2 = MultilayerGraph(W_2, connect = Multiplex())
	partition_2 = SpectralPartition(mlgraph_2)
	seba_part_2 = SEBAPartition(partition_2,3)
end

# ╔═╡ 158bf43f-cc38-46d6-980d-6cfe75522c4c
# Plot the network at its first, middle and last time steps.

# ╔═╡ 0dbbc497-606f-4ba6-9d7c-fbfff6c73853
begin
	p = [plot(mlgraph_2,x) for x in [1,32,60]]
	plot(p...)
end

# ╔═╡ f6ece308-8af8-4926-a4c3-a67c6da63456
# Plot leading eigenvectors of the supra-Laplacian

# ╔═╡ 3c421c9d-28d0-48c2-983f-225982a5be08
plot(plot(partition_2)[1][1:6]...)

# ╔═╡ b0640fd6-7d34-4336-b585-6042650172ba
# Plot SEBA vectors (SEBA performed with doubling) corresponding to the first three  spatial evecs of the supra-Laplacian.

# ╔═╡ 41e7fdb1-e74b-4969-a582-551b2c4e63fc
plot(plot(seba_part_2)[1]...)

# ╔═╡ Cell order:
# ╠═980691f0-3ec3-11ef-0ef1-e7653399f7de
# ╠═7324ec91-a825-47c4-b2aa-2d8a22b5ba3c
# ╠═b8449181-55d3-42e1-83cf-57ede95e1033
# ╠═3f621722-cced-45f8-81ce-485234c289a3
# ╠═8f50edf1-fe0c-4e43-8fa6-d02639c22416
# ╠═dd9c3a1e-6bad-4895-8342-33b6fa9a43aa
# ╠═cb73e5c2-b133-430e-9ca3-d052b7262bba
# ╠═15b20c67-a92b-401c-81ff-11e7209afdd2
# ╠═daa0159f-da2b-4bb3-8abc-51ed89124b25
# ╠═9d67cfd0-d114-4bc9-b765-afe8b2221469
# ╠═5693c8ae-c7cb-425f-91b0-5ac8d38a1c55
# ╠═158bf43f-cc38-46d6-980d-6cfe75522c4c
# ╠═0dbbc497-606f-4ba6-9d7c-fbfff6c73853
# ╠═f6ece308-8af8-4926-a4c3-a67c6da63456
# ╠═3c421c9d-28d0-48c2-983f-225982a5be08
# ╠═b0640fd6-7d34-4336-b585-6042650172ba
# ╠═41e7fdb1-e74b-4969-a582-551b2c4e63fc
