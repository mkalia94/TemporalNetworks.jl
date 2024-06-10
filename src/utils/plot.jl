import Plots: plot

function plot(mlgraph :: MultilayerGraph{T}, layer :: Int64; labels = nothing, kwargs...) where T <: Union{Multiplex, NonMultiplex}

    ind_active = T <: NonMultiplex ? find_active(mlgraph)[layer] : 1:mlgraph.N
    W(x) = mlgraph.W[x][ind_active, ind_active]
    nm = mlgraph.N < 10 ? ["$(i) " for i in 1:mlgraph.N][ind_active] : [["$(i) " for i in 1:9]; ["$(i)" for i in 10:mlgraph.N]][ind_active]

    edge_colors = ones(size(W(layer))) .|> Int64
    edge_color!(W, mlgraph.T, layer, edge_colors) 
    node_colors = isnothing(labels) ? Int64.(ones(length(ind_active))) :  color_vector(Int64.(labels[:,layer]))[ind_active]
    
    g = graphplot(W(layer), names=nm,nodesize=0.2,method = :circular, nodeshape = :circle,  edgecolor = edge_colors,  nodecolor=node_colors, curvature = 0.01, kwargs...)
    plot(g, dpi=300)
end 


function edges_to_delete(W, layer :: Int64, edge_colors :: Matrix{Int64})
    edges = [z for z ∈ edgelist(W(layer)) if z ∉ edgelist(W(layer+1)) ]
    map(x->edge_colors[x...]=2, edges)
end

function edges_to_delete(W :: Matrix{Float64}, layer :: Int64)
    edges = [z for z ∈ edgelist(W[layer]) if z ∉ edgelist(W[layer+1]) ]
end

function edges_added(W, layer :: Int64, edge_colors :: Matrix{Int64})
    edges = [z for z ∈ edgelist(W(layer)) if z ∉ edgelist(W(layer-1)) ]
    map(x->edge_colors[x...]=3, edges)
end

function edges_added(W ::  Matrix{Float64}, layer :: Int64)
    edges = [z for z ∈ edgelist(W[layer]) if z ∉ edgelist(W[layer-1]) ]
end

function edge_color!(W , T :: Int64, layer :: Int64, edge_colors :: Matrix{Int64})
    layer == 1 && return edges_to_delete(W, layer, edge_colors)
    layer == T && return edges_added(W, layer, edge_colors)
    edges_to_delete(W, layer, edge_colors)
    edges_added(W, layer, edge_colors)
end

function plot(partition :: SpectralPartition{MultilayerGraph{T}, N}; vecs = nothing, clims = nothing, kwargs...) where {T,N}

    h = []
    vec_inds, vecs_to_plot = _vecs_to_plot(partition, vecs = vecs)
    
    for (i, ind) in enumerate(vec_inds)
        yy = maximum(abs.(partition.evecs[:,ind]))
        _clims = isnothing(clims) ?  (-yy, yy) : nothing
        push!(h,heatmap(reshape(vecs_to_plot[:,i], partition.graph.N, partition.graph.T),title="Orig. Evec $(ind)",c=cgrad(:RdBu), dpi=300, clims = (-yy,yy), kwargs...))
    end
    h
end

function _vecs_to_plot(partition :: SpectralPartition{MultilayerGraph{T}, N}; vecs = nothing, clims = nothing, kwargs...) where {T <: NonMultiplex, N <: Normalization}
    vecs_to_plot = isnothing(vecs) ? Array(1:10) : vecs
    return vecs_to_plot, partition.evecs[:, vecs_to_plot]
end

function _vecs_to_plot(partition :: SpectralPartition{MultilayerGraph{T}, N}; vecs = nothing, clims = nothing, kwargs...) where {T <: Multiplex, N <: Normalization}
    vecs_to_plot = isnothing(vecs) ? Array(1:10) : vecs
    return vecs_to_plot, partition.evecs[:, vecs_to_plot]
end


function plot(spart:: SEBAPartition{MultilayerGraph{T}, N}; kwargs...) where {T, N}
    h = []
    inds = 1:size(spart.vecs)[2]
    for i in inds
        push!(h,heatmap(reshape( spart.vecs[:,i],spart.partition.graph.N,spart.partition.graph.T),title="SEBA $(i), cut: $(round(spart.cuts[i], digits=3))",c=:tempo, xticks=get_ticks(spart.partition.graph.T), clims=(0,1), kwargs...))
    end
    h
end


