using Pkg; Pkg.activate("../.")
using TemporalNetworks, Plots, FileIO, JLD2

W = FileIO.load("senator.jld2","W")
labels_party = FileIO.load("senator.jld2","labels_party")
labels_state = FileIO.load("senator.jld2","labels_state")

mlgraph= MultilayerGraph(W, connect = NonMultiplexCompressed())
partition = SpectralPartition(mlgraph, compute_a = RayleighBalancing(2))

# Reordering vertices for clarity
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

# Plots, with state results
W_state = FileIO.load("state.jld2", "W")
mlgraph_state = MultilayerGraph(W_state)
partition_state = SpectralPartition(mlgraph_state)


