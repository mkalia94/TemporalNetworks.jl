### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 812e0004-32d8-11ef-37e6-1d4b240c9a4b
begin
	using Pkg; Pkg.activate("../.")
	using TemporalNetworks, Plots, FileIO, JLD2
end

# ╔═╡ c558534e-c1e5-41af-93b2-69a16ebe0af9
# Load packages for spectral partitioning, plotting and loading datasets.

# ╔═╡ 4a7d678b-ffc0-4ece-8bbc-2eadcd02926d
# Load datasets and construct graph and spectral partition instances.

# ╔═╡ d96268dd-5dcc-4a10-8f89-fbe12e9765c9
begin
	W = FileIO.load("senator.jld2","W")
	labels_party = FileIO.load("senator.jld2","labels_party")
    labels_state = FileIO.load("senator.jld2","labels_state")
	
	mlgraph= MultilayerGraph(W, connect = NonMultiplexCompressed())
	partition = SpectralPartition(mlgraph, compute_a = RayleighBalancing(2))
end

# ╔═╡ 65b1a46e-19ca-458d-9bbe-6b52982248af
begin
	using LinearAlgebra
	idl = partition.L_spat + partition.a^2 .* partition.L_temp
    min_deg = minimum(diag(idl))
    max_deg = maximum(diag(idl))
	ind_new = findall(x->(x ≥ max_deg || x ≤ min_deg), partition.evals)
end


# ╔═╡ 9ac21e96-8cb4-4f73-8605-637247ed99ef
# Now we plot eigenvectors of the inflated dynamic Laplacian to select vectors which will induce the desired spacetime partition. 

# ╔═╡ 5a598148-8abe-41ef-9239-2465a1b66231
h = plot(partition); plot(h[1][1:4]..., dpi=300)

# ╔═╡ a791feb8-88d4-47d4-81ea-c9b8faaf8497
# Clearly evec 3 is the first `spatial' eigenvector, which we use to construct the spacetime Laplacian. First, we carefully reorder vertices according to their temporal mean values, for visual convenience.

# ╔═╡ d4830a3f-cf09-4bf4-8967-3affcbafcdb4
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

# ╔═╡ 4895258b-315f-4d3f-8d48-a8cf228b7b3e
# We now use v and ordering from the previous cell to reorder eigenvectors, labels and the network adjacencies appropriately.

# ╔═╡ 2bc2b626-f650-40b3-8252-1d979e03d5ac
begin
	partition.evecs = partition.evecs[ordering, :]
	lab_party = labels_party[v]
	lab_state = labels_state[v]
    mlgraph.W = [x[v,v] for x in mlgraph.W]
end

# ╔═╡ 3761f001-40a7-4b13-ab4e-3f5bfb2a20d1
# We do a small preprocessing step for the party affiliations

# ╔═╡ 424a0bcf-75aa-424b-a8cc-6ee1ee693313
	for i in 1:mlgraph.N
    	if lab_party[i] == 328
        	lab_party[i] = 0
    	elseif lab_party[i] == 100
        	lab_party[i] = 1
    	elseif lab_party[i] == 200
        	lab_party[i] = -1
    	end
	end

# ╔═╡ 82b79d39-8910-4f1d-9c04-cdbfbe63a060
# Now we load the state network data and create its spectral partitioning instance

# ╔═╡ e531470d-98e3-48f8-a086-35ec654a5046
begin
	    W_state = FileIO.load("state.jld2", "W");
	    lab_statelist = FileIO.load("state.jld2", "labels");
		mlgraph_state = MultilayerGraph(W_state);
		partition_state = SpectralPartition(mlgraph_state)
end

# ╔═╡ 77898ce1-0e53-429d-822c-d34af5462ca9
# Next, similar to the senator case, we plot eigenvectors of the inflated dynamic Laplacian corresponding to the state network.

# ╔═╡ 14d35b45-a8a0-49d0-9c0d-4eeebbf04a45
h_state = plot(partition_state); plot(h_state[1][1:4]...)

# ╔═╡ bf0d84f3-8246-412c-a4df-bd3e63e6c597
# Clearly the first spatial eigenvector is the second one (note that this is a multiplex network as opposed to the senator case). Thus we use the second eigenvector for spectral partitioning. We reorder the state vertices as well, by the magnitude of the eigenvector at the last time step.

# ╔═╡ 5e0b55ab-9271-4649-b0ca-cf0a602fe98e
begin
	evec_crit_state = reshape(partition_state.evecs[:,2],partition_state.graph.N, partition_state.graph.T)
	v_state = sortperm(evec_crit_state[:,end])
end

# ╔═╡ e433a557-667d-452b-b4fd-9d8f8c456dec
# Finally, in order to compare the senator and state results, we compute the projection of the third eigenvector of the senator network onto the state network. This is done by finding the two senators corresponding to each state per time step, and assigning them the value obtained from the senator network eigenvector.

# ╔═╡ ca5b9847-613f-4b2d-907d-9c1862bdf5d7
begin
	evecs_temp = [reshape(partition.evecs[:,3], mlgraph.N, mlgraph.T),]
    sens_ = find_active(mlgraph) 
	state_proj = [zeros(2*length(lab_statelist[v_state]), partition.graph.T) for i in 1:length(evecs_temp)]
	for i in 1:length(state_proj)
    	for x in 1:length(lab_statelist[v_state]), y in 1:size(state_proj[i])[2]
       	 	ind_ = findall(z->(z ∈ sens_[y] && lab_state[z] == lab_statelist[v_state][x]), 1:partition.graph.N)
       	 	if length(ind_)==1
        	    state_proj[i][(x-1)*2+1,y] = evecs_temp[i][ind_[1],y]
       		else
        	    state_proj[i][(x-1)*2+1:2*x,y] = evecs_temp[i][ind_[1:2], y]
        	    # None of this accounts for states that may have more than two senators in a single congress (with incomplete terms)
       	    end
    	end
	end
	# Filter
	for i in 1:length(state_proj)
    	for (j,x) in enumerate(state_proj[i])
        	if x == 0.0
            	state_proj[i][j] = NaN
        	end
    	end
	end
	state_proj
end

# ╔═╡ 05e5f466-2c63-4bb2-af97-ace816b72a32
# The object `state_proj' can now be compared with the second eigenvector of the state_network

# ╔═╡ 3b358fd4-15a5-401c-9fed-9128aceefc48
# (First, an appropriate scaling for plotting)

# ╔═╡ a5b645f8-0525-4752-83d2-60818f05f852
begin
	state_proj[1] .*= sqrt(200);
	partition_state.evecs .*= sqrt(200);
end

# ╔═╡ d7dd157a-97a9-4821-8361-92defca202dc
begin	
	labelss = vcat([[x; ""] for x in lab_statelist[v_state]]...);
	
	h1 = heatmap(reshape(partition_state.evecs[:,2],mlgraph_state.N,mlgraph_state.T)[v_state,:], c=cgrad(:RdBu), size=(600,800), yticks=(1:50,lab_statelist[v_state]), xticks=(1:2:11,1987:4:2009), dpi=300, clims = (-1,1))
	
	h2 = heatmap(-state_proj[1], c=cgrad(:RdBu), size=(600,800), yticks=(1:100,labelss), xticks=(1:2:11,1987:4:2009), dpi=300, grid=false, clims=(-1,1))

	hline!(h2,0.5:2.0:(101+0.5), c=:black, labels="", lw=1)
	
	plot(h1, h2, size=(900,700), dpi=300)
	savefig("proj.png"); gui()
	plot(h1, h2, size=(900,700), dpi=300)
end



# ╔═╡ 8a088ae0-c0e3-4583-baca-8b08a41c25a1
[isspat(mlgraph, partition.norm, partition.L_spat, partition.L_temp, partition.a, x, nothing) for x in 1:15]

# ╔═╡ cdb839be-a89a-4c55-91ba-7f712f69ad7d
s1 = SEBAPartition(partition, [3])

# ╔═╡ f5156ea8-6815-42e6-9ec1-511da0b0c845
s1.cuts

# ╔═╡ Cell order:
# ╠═c558534e-c1e5-41af-93b2-69a16ebe0af9
# ╠═812e0004-32d8-11ef-37e6-1d4b240c9a4b
# ╠═4a7d678b-ffc0-4ece-8bbc-2eadcd02926d
# ╠═d96268dd-5dcc-4a10-8f89-fbe12e9765c9
# ╠═9ac21e96-8cb4-4f73-8605-637247ed99ef
# ╠═5a598148-8abe-41ef-9239-2465a1b66231
# ╠═a791feb8-88d4-47d4-81ea-c9b8faaf8497
# ╠═d4830a3f-cf09-4bf4-8967-3affcbafcdb4
# ╠═4895258b-315f-4d3f-8d48-a8cf228b7b3e
# ╠═2bc2b626-f650-40b3-8252-1d979e03d5ac
# ╠═3761f001-40a7-4b13-ab4e-3f5bfb2a20d1
# ╠═424a0bcf-75aa-424b-a8cc-6ee1ee693313
# ╠═82b79d39-8910-4f1d-9c04-cdbfbe63a060
# ╠═e531470d-98e3-48f8-a086-35ec654a5046
# ╠═77898ce1-0e53-429d-822c-d34af5462ca9
# ╠═14d35b45-a8a0-49d0-9c0d-4eeebbf04a45
# ╠═bf0d84f3-8246-412c-a4df-bd3e63e6c597
# ╠═5e0b55ab-9271-4649-b0ca-cf0a602fe98e
# ╠═e433a557-667d-452b-b4fd-9d8f8c456dec
# ╠═ca5b9847-613f-4b2d-907d-9c1862bdf5d7
# ╠═05e5f466-2c63-4bb2-af97-ace816b72a32
# ╠═3b358fd4-15a5-401c-9fed-9128aceefc48
# ╠═a5b645f8-0525-4752-83d2-60818f05f852
# ╠═d7dd157a-97a9-4821-8361-92defca202dc
# ╠═65b1a46e-19ca-458d-9bbe-6b52982248af
# ╠═8a088ae0-c0e3-4583-baca-8b08a41c25a1
# ╠═cdb839be-a89a-4c55-91ba-7f712f69ad7d
# ╠═f5156ea8-6815-42e6-9ec1-511da0b0c845
