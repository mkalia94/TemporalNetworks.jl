using Pkg
Pkg.activate(".")

using TemporalNetworks, Plots

list = [0,1]
clusters = [[Array(1:5)],
            [Array(1:3), Array(4:5)]]
degrees = [[4],[2,1]]
η = 0.8
N = 5
T = 10

block = BlockGraph(N, T, list, η, clusters, degrees)
W1 = block()

vals = range(0.3,5.0,length=200)
for (i,W) in enumerate(W1)
    for (j,x) = enumerate(W)
        W1[i][j] = x*rand(vals)
    end
end

W1 = [(x+x')/2 for x in W1]

mlgraph = MultilayerGraph(W1)
@info "T = $(mlgraph.T)"
partition = SpectralPartition(mlgraph)
partition_norm = SpectralPartition(mlgraph, partition.a, norm = DegreeNormalization())
p1 = plot(partition.evals[1:20] ./ partition.a^2 ./2)
plot!(partition_norm.evals[1:20])


