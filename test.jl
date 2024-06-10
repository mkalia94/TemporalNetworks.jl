using Pkg; Pkg.activate(".")
using TemporalNetworks
# Block graphs
list = [0,1,2]
degrees = nothing
η = 0.8
clusters = nothing
block = BlockGraph(20, 10, list, η, clusters, degrees)
W1 = block()

# Block graphs nonmultiplex
η = 0.8
list = [2,1]
clusters = [[Array(1:5), Array(6:20), Array(21:25)],
                [Array(1:10), Array(11:20)]] # you can ignore this line
degrees = [[4,6,4],
        [6,6]] # Ignore this line
evolve = 1
block = BlockGraphNonMultiplex(25, 10, list, η, clusters, degrees, evolve)
W2 = block() |> Vector{Matrix{Float64}}


# Tests for non-multiplex networks
mlgraph_nonmultiplex = MultilayerGraph(W2, connect = NonMultiplexCompressed())
partition_nonmultiplex = SpectralPartition(mlgraph_nonmultiplex, compute_a = RayleighBalancing(3)) # rayleigh balancing on third eigenvalue
# partition_nonmultiplex2 = SpectralPartition(mlgraph_nonmultiplex, compute_a = SpatTempMatchingNonMultiplex(5, nothing, 0.2)) # match evec 3, pick best temporal eigenvector, threshold for inner products is 0.1 
# The above code is commented out as it somehow doesn't work for the nonmultiplex block graph
seba_part_nonmultiplex = SEBAPartition(partition_nonmultiplex)

# Tests for multiplex networks
mlgraph = MultilayerGraph(W1, connect = Multiplex())
partition_1 = SpectralPartition(mlgraph) # Classic eigenvalue matching 
partition_2 = SpectralPartition(mlgraph, compute_a = RayleighBalancing(2)) # Rayleigh Balancing on second eigenvalue
seba_part  = SEBAPartition(partition_1)

p1 = plot(mlgraph, 1)
p2 = plot(mlgraph_nonmultiplex, 1)
p3 = plot(partition_nonmultiplex)
p4 = plot(partition_1)
p5 = plot(seba_part)
p6 = plot(seba_part_nonmultiplex)

