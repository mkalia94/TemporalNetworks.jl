using Pkg
Pkg.activate(".")

using TemporalNetworks, Plots, LinearAlgebra

list = [0,1,0]
clusters = [[Array(1:20)],
            [Array(1:10), Array(11:20)],
            [Array(1:20)]]
degrees = [[14],[6,6],[14]]
η = 0.8
N = 20
T = 20

block = BlockGraph(N, T, list, η, clusters, degrees)
W1 = block()


struct RandomMultiplex <: TemporalNetworks.TemporalConnectivity end

function (ff :: RandomMultiplex)(mlgraph :: MultilayerGraph)
    Wt = hcat([[rand(range(0,2,length=200)) for i in 1:mlgraph.T] for i in 1:mlgraph.T]...) 
    Wt = (Wt + Wt')/2
    TemporalNetworks.directsum(mlgraph.W),  kron(Wt, Matrix{Float64}(I, mlgraph.N, mlgraph.N))
end


function _plot_eigs(partition :: SpectralPartition{MultilayerGraph{T}, N}, R) where {T <: RandomMultiplex, N}
    nothing
end

mlgraph = MultilayerGraph(W1, connect = RandomMultiplex())

partition = SpectralPartition(mlgraph, compute_a = RayleighBalancing(2))
partition_norm = SpectralPartition(mlgraph, partition.a, norm = DegreeNormalization())
p1 = plot(partition.evals[1:20] ./ partition.a^2 ./2)
plot!(partition_norm.evals[1:20])

