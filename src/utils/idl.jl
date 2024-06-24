function DynamicLaplacian(mlgraph :: MultilayerGraph)
    norm(sum(mlgraph.W) ./ mlgraph.T)
end

function DynamicLaplacian(sp :: SpectralPartition)
    DynamicLaplacian(sp.graph)
end

function InflatedDynamicLaplacian(sp :: SpectralPartition)
    norm(sp.L_spat + sp.a^2 .* sp.L_temp)
end
