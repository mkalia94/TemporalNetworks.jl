using Pkg; Pkg.activate(".")
using TemporalNetworks


#Notes: using MultilayerGraph, SpectralPartition, SEBAPartition we 
# construct class instances (like in C) that contain cool info
# 1. For MultilayerGraph instances, say ml, ml.graph gives the graph, ml.connect gives a function that connects the list of adjacencies, this can be nonmultiplex compressed or multiplex
# 2. For SpectralPartition instances, say sp, sp.norm gives how the Laplacian is normalized, sp.evecs gives eigenfunctions (lifted if nonmultiplex), sp.evals gives eigenvalues etc.
# 3. For SEBAPartition instances, say seba, seba.inds gives the eigenvector indices used to create the SEBA partition, seba.vecs gives the SEBA vectors and seba.cuts gives the corresponding cut values based on normalization. 
# All objects can be plotted using plot(...) and they will plot the corresponding graph objects. Plotting ml returns the graph, plotting sp returns the evecs and plotting seba returns the SEBA vectors with cut values.  

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


