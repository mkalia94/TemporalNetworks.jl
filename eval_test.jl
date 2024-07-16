using Pkg
Pkg.activate(".")

using TemporalNetworks, Plots, LinearAlgebra

# Parameters for block graphs
list = [0,1,2]
degrees = nothing
η = 0.8
clusters = nothing

N = 10:20
T = 20:30

counts = []

for i in 1:100
    @show i
    block = BlockGraph(rand(N), rand(T), list, η, clusters, degrees)
    W1 = block()
    mlgraph = MultilayerGraph(W1)
    try
        partition = SpectralPartition(mlgraph)
        idl = partition.L_spat + partition.a^2 .* partition.L_temp
        min_deg = minimum(diag(idl))
        max_deg = maximum(diag(idl))
        @show ind = findall(x-> (x <= min_deg), partition.evals)
        push!(counts, length(ind)/mlgraph.N/mlgraph.T)
    catch e
        @show e
        continue
    end
end

    
