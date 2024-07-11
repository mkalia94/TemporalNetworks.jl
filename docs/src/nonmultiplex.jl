using Pkg; Pkg.activate("../../.")
using TemporalNetworks, Plots, LaTeXStrings

η = 0.8
list = [2,1]
clusters = [[Array(1:5), Array(6:20), Array(21:25)],
                [Array(1:10), Array(11:20)]] # you can ignore this line
degrees = [[4,6,4],
        [6,6]] # Ignore this line
evolve = 1
block = BlockGraphNonMultiplex(25, 10, list, η, clusters, degrees, evolve)
W2 = block() |> Vector{Matrix{Float64}}


mlgraph_nonmultiplex = MultilayerGraph(W2, connect = NonMultiplexCompressed())
partition_nonmultiplex = SpectralPartition(mlgraph_nonmultiplex, compute_a = RayleighBalancing(3)) # Rayleigh balancing on third eigenvalue
seba_part_nonmultiplex = SEBAPartition(partition_nonmultiplex,[2,])

plot(plot(partition_nonmultiplex)[1][1:6]..., size = (700,400), dpi=300)
savefig("figs/evecs_nonm.png")

plot(plot(partition_nonmultiplex)[2], size = (350,250), dpi=300)
savefig("figs/evals_nonm.png")

plot(plot(seba_part_nonmultiplex)..., size=(700,400), dpi=300)
savefig("figs/SEBA_nonm.png")

plot(seba_part_nonmultiplex,[1,2], dpi=300)
savefig("figs/heatmap_nonm.png")

xx = leiden_slice(partition_nonmultiplex)
heatmap(xx[2], c=cgrad([:white, :orange, :red]), size=(400,300), dpi=300)
savefig("figs/leiden-slice-nonm.png")
