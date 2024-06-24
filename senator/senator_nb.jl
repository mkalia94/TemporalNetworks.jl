### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ 3ec29f89-f8e4-4b02-8151-2b19f9694785
using Pkg; Pkg.activate("../.")

# ╔═╡ 00c79aac-1670-484c-b516-0cd2ed2da044
# Load senator network files and construct inflated dynamic Laplacian

# ╔═╡ 06468a08-bca2-4ab4-b8b4-e789da5dfe51
h = plot(partition); plot(h[1:5]..., size=(900,500))

# ╔═╡ d308880d-bfe7-4501-a14c-dda42711c31d
# From the above plot, it is clear that the first `spatial' eigenvector is the third one. We use this eigenvector to carefully reorder the vertices for clarity.

# ╔═╡ b9d18a02-722b-46e0-bf83-123eaab88cd1
begin
	evec_crit = reshape(partition.evecs[:,3],mlgraph.N, mlgraph.T)
	sens = find_active(mlgraph)

	lengths = zeros(mlgraph.N)
	for i in 1:mlgraph.N
    	ctr = 0
    	for x in sens
        	ind = findall(y->y==i,x)
        	ctr += length(ind)
    	end
    	lengths[i] = ctr
	end

	aa = sum(hcat([evec_crit[:,i] for i in 1:mlgraph.T]...), dims=2) ./ lengths
	v = sortperm(aa, dims=1)[:,1]
	ordering = vcat([v .+ (i-1)*mlgraph.N for i in 1:mlgraph.T]...)
end

# ╔═╡ bd4ef6e3-9aba-4c84-ba73-10eda482bb98
# We rearrange the vertices and the corresponding labels. For the vertices we do so only in the spectral partition instance.

# ╔═╡ cddd6ebb-601e-4d02-bc32-aa85269d16c3
# Now we load the network of states and compute its spectral partition

# ╔═╡ 1f207d66-5cc2-4272-86bc-4ef5007f2071
begin
	W_state = FileIO.load("states.jld2", "W");
	labels_statelist = FileIO.load("states.jld2", "labels");
	mlgraph_state = MultilayerGraph(W_state);
	partition_state = SpectralPartition(mlgraph_state)
end

# ╔═╡ 1cc0afcd-a1e4-47a1-9158-a3c2c2a38da1
h_state = plot(partition_state); plot(h_state[1:4]...)

# ╔═╡ c178d5f8-33e4-4ecf-8609-d0110cfebf45
# From the above plot, it is clear that evec 2 is the first spatial eigenvector. We rearrange the state vertices such that the projections reflect a trend of strong democrat to strong republican regimes.

# ╔═╡ 9df8fd1e-c239-4116-8a91-083e743a8df9
begin
	evec_crit_state = reshape(partition_state.evecs[:,2],partition_state.graph.N, partition_state.graph.T)
	v_state = sortperm(evec_crit_state[:,end])
end

# ╔═╡ a59e9b25-9b02-4844-836c-12385fe7a6f0
begin
	include("project_sens_to_state.jl")
	state_proj = project_sens_to_state(partition, [reshape(partition.evecs[:,3], mlgraph.N, mlgraph.T),], labels_state, labels_statelist[2:end][v_state])
	state_proj[1] .*= sqrt(200)
	partition_state.evecs .*= sqrt(200)
end

# ╔═╡ a45de2ee-7cd3-4404-96be-a19279d9cb5b
# Now we project the senator results onto the state using the project_sens_to_state function supplied with the notebook. For each state, there exist minimum two senators per time step. Thus we obtain a multiplex network of 100 senators per time step, to be compared with the state network.

# ╔═╡ 4e5da72e-a482-47f9-b1e0-fc81781bb57f
# Finally we are ready to plot senator and state results side-by-side

# ╔═╡ 34906585-dc9a-4b70-982e-2c79eed68d82
begin
	labelss = vcat([[x; ""] for x in labels_statelist[2:end][v_state]]...)
	h1 = heatmap(reshape(-partition_state.evecs[:,2],mlgraph_state.N,mlgraph_state.T)[v_state,:], c=cgrad(:RdBu), size=(600,800), yticks=(1:50,labels_statelist[2:end][v_state]), xticks=(1:2:11,1987:4:2009), dpi=300, clims = (-0.84,0.84))
	h2 = heatmap(state_proj[1], c=cgrad(:RdBu), size=(600,800), yticks=(1:100,labelss), xticks=(1:2:11,1987:4:2009), dpi=300, grid=false, clims=(-1,1))
	hline!(h2,0.5:2.0:(101+0.5), c=:black, labels="", lw=1)
	plot(h1, h2, size=(900,700))
	
end

# ╔═╡ cf22331e-7cf3-4e9b-be85-91eb80553692
begin
	using TemporalNetworks, Plots, FileIO, JLD2
	
	W = FileIO.load("senator.jld2","W")
	labels_party = FileIO.load("senator.jld2","labels_party")
	labels_state = FileIO.load("senator.jld2","labels_state")
	
	mlgraph= MultilayerGraph(W, connect = NonMultiplexCompressed())
	partition = SpectralPartition(mlgraph, compute_a = RayleighBalancing(2))
end

# ╔═╡ 35021c84-8c90-4c1b-941d-3c1b5e4d008e
begin
	partition.evecs = partition.evecs[ordering, :]
	labels_party = labels_party[v]
	labels_state = labels_state[v]


	for i in 1:mlgraph.N
    	if labels_party[i] == 328
        	labels_party[i] = 0
    	elseif labels_party[i] == 100
        	labels_party[i] = 1
    	elseif labels_party[i] == 200
        	labels_party[i] = -1
    	end
	end


	labels_plot = deepcopy(evec_crit[v,:])

	for i in 1:mlgraph.N, j in 1:mlgraph.T
    	if labels_plot[i,j] != 0
        	labels_plot[i,j] = labels_party[i]
    	end
	end
end

# ╔═╡ Cell order:
# ╠═00c79aac-1670-484c-b516-0cd2ed2da044
# ╠═3ec29f89-f8e4-4b02-8151-2b19f9694785
# ╠═cf22331e-7cf3-4e9b-be85-91eb80553692
# ╠═06468a08-bca2-4ab4-b8b4-e789da5dfe51
# ╠═d308880d-bfe7-4501-a14c-dda42711c31d
# ╠═b9d18a02-722b-46e0-bf83-123eaab88cd1
# ╠═bd4ef6e3-9aba-4c84-ba73-10eda482bb98
# ╠═35021c84-8c90-4c1b-941d-3c1b5e4d008e
# ╠═cddd6ebb-601e-4d02-bc32-aa85269d16c3
# ╠═1f207d66-5cc2-4272-86bc-4ef5007f2071
# ╠═1cc0afcd-a1e4-47a1-9158-a3c2c2a38da1
# ╠═c178d5f8-33e4-4ecf-8609-d0110cfebf45
# ╠═9df8fd1e-c239-4116-8a91-083e743a8df9
# ╠═a45de2ee-7cd3-4404-96be-a19279d9cb5b
# ╠═a59e9b25-9b02-4844-836c-12385fe7a6f0
# ╠═4e5da72e-a482-47f9-b1e0-fc81781bb57f
# ╠═34906585-dc9a-4b70-982e-2c79eed68d82
