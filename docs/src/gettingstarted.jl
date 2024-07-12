using Pkg; Pkg.activate("../../.")
using TemporalNetworks, Plots, LaTeXStrings
list = [0,2]; degrees = nothing; η = 0.8; 
clusters = [[Array(1:20)],[Array(1:10), Array(11:15), Array(16:20)] ]
degrees = [[14],[5,4,4]]
block = BlockGraph(20, 15, list, η, clusters, degrees)
W1 = block();
mlgraph = MultilayerGraph(W1, connect = Multiplex())
partition = SpectralPartition(mlgraph)
seba_part = SEBAPartition(partition,2)

p = [plot(mlgraph,x) for x in 1:5:mlgraph.T]
plot(p..., dpi=300, size = (700,400))
savefig("figs/graph.svg")

plot(plot(partition)[1][1:6]..., size = (700,400), dpi=300)
savefig("figs/evecs.svg")

plot(plot(partition)[2], size = (350,250), dpi=300)
savefig("figs/evals.svg")

plot(plot(seba_part)[2], size=(400,300), dpi=300)
savefig("figs/Ratios.svg")

plot(plot(seba_part)[1]..., size=(700,400), dpi=300)
savefig("figs/SEBA.svg")

plot(seba_part,[1,2,3], dpi=300)
savefig("figs/heatmap.svg")

xx = leiden_slice(partition)
heatmap(xx[2], c=cgrad([:white, :orange, :red]), size=(400,300), dpi=300)
savefig("figs/leiden-slice.svg")
