import Plots: plot

function plot(mlgraph :: MultilayerGraph{T}, layer :: Int64; labels = nothing, kwargs...) where T <: Union{Multiplex, NonMultiplex}

    ind_active = T <: NonMultiplex ? find_active(mlgraph)[layer] : 1:mlgraph.N
    W(x) = mlgraph.W[x][ind_active, ind_active]
    nm = mlgraph.N < 10 ? ["$(i) " for i in 1:mlgraph.N][ind_active] : [["$(i) " for i in 1:9]; ["$(i)" for i in 10:mlgraph.N]][ind_active]

    edge_colors = ones(size(W(layer))) .|> Int64
    edge_color!(W, mlgraph.T, layer, edge_colors) 
    node_colors = isnothing(labels) ? Int64.(ones(length(ind_active))) :  color_vector(Int64.(labels[:,layer]))[ind_active]
    
    g = graphplot(W(layer), names=nm,nodesize=0.2,method = :circular, nodeshape = :circle,  edgecolor = edge_colors,  nodecolor=node_colors, curvature = 0.01, title="t=$(layer)", kwargs...)
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

function plot(partition :: SpectralPartition{MultilayerGraph{T}, N}, R :: Int64; vecs = nothing, clims = nothing, kwargs...) where {T,N}
    _plot_partition(partition, vecs = vecs, clims =clims, kwargs...), _plot_eigs(partition, R)
end

function plot(partition :: SpectralPartition{MultilayerGraph{T}, N}; vecs = nothing, clims = nothing, kwargs...) where {T,N}
    _plot_partition(partition, vecs = vecs, clims =clims, kwargs...), _plot_eigs(partition, 20)
end

function _plot_partition(partition :: SpectralPartition{MultilayerGraph{T}, N}; vecs = nothing, clims = nothing, kwargs...) where {T,N}

    h = []
    vec_inds, vecs_to_plot = _vecs_to_plot(partition, vecs = vecs)
    
    for (i, ind) in enumerate(vec_inds)
        yy = maximum(abs.(partition.evecs[:,ind]))
        _clims = isnothing(clims) ?  (-yy, yy) : nothing
        push!(h,heatmap(reshape(vecs_to_plot[:,i], partition.graph.N, partition.graph.T),title="Orig. Evec $(ind)",c=cgrad(:RdBu), dpi=300, clims = (-yy,yy), kwargs...))
    end
    h
end

function _plot_eigs(partition :: SpectralPartition{MultilayerGraph{T}, N}, R) where {T <: Multiplex, N}
    labels_eigs = zeros(R)
    ind = findall(x->isspat(partition.evecs[:,x],partition.graph.N,partition.graph.T,0.01) == true, 1:R)
    ind_not = [x for x in 1:R if x ∉ ind]
    p5 = scatter(ind, partition.evals[ind], color=:red, label=L"\Lambda_{k,a}^{\rm spat}", xlabel=L"k", title=L"\textrm{Eigenvalues } \ \Lambda_{k,a} \ \textrm{ vs } \ k ", titlefontsize=10)
    scatter!(p5, ind_not, partition.evals[ind_not], color=:black, markershape=:square, label=L"\Lambda_{k,a}^{\rm temp}", dpi=300, legend=:bottomright)
    p5
end

function _plot_eigs(partition :: SpectralPartition{MultilayerGraph{T}, N}, R) where {T <: NonMultiplex, N}
    labels_eigs = zeros(R)
    ind = findall(x->isspat(partition.graph, partition.norm, partition.L_spat, partition.L_temp, partition.a, x, nothing) < 0.01, 1:R)
    ind_not = [x for x in 1:R if x ∉ ind]
    p5 = scatter(ind, partition.evals[ind], color=:red, label=L"\Lambda_{k,a}^{\rm spat}", xlabel=L"k", title=L"\textrm{Eigenvalues } \ \Lambda_{k,a} \ \textrm{ vs } \ k ", titlefontsize=10)
    scatter!(p5, ind_not, partition.evals[ind_not], color=:black, markershape=:square, label=L"\Lambda_{k,a}^{\rm temp}", dpi=300, legend=:bottomright)
    p5
end

function _vecs_to_plot(partition :: SpectralPartition{MultilayerGraph{T}, N}; vecs = nothing, clims = nothing, kwargs...) where {T <: NonMultiplex, N <: Normalization}
    vecs_to_plot = isnothing(vecs) ? Array(1:10) : vecs
    return vecs_to_plot, partition.evecs[:, vecs_to_plot]
end

function _vecs_to_plot(partition :: SpectralPartition{MultilayerGraph{T}, N}; vecs = nothing, clims = nothing, kwargs...) where {T <: Multiplex, N <: Normalization}
    vecs_to_plot = isnothing(vecs) ? Array(1:10) : vecs
    return vecs_to_plot, partition.evecs[:, vecs_to_plot]
end


function plot(spart :: SEBAPartition{MultilayerGraph{T}, N}; kwargs...) where {T<: Multiplex, N}
    _plot_SEBA(spart), _plot_SEBA_trend(spart, 5)
end

function plot(spart :: SEBAPartition{MultilayerGraph{T}, N}, R :: Int64; kwargs...) where {T <: Multiplex, N}
    _plot_SEBA(spart), _plot_SEBA_trend(spart, R)
end

function plot(spart :: SEBAPartition{MultilayerGraph{T}, N}; kwargs...) where {T<: NonMultiplex, N}
    _plot_SEBA(spart)
end

function plot(spart :: SEBAPartition{MultilayerGraph{T}, N}, good_SEBA :: Vector{Int64}; kwargs...) where {T, N}
    p, _, _ = _plot_SEBA_heatmapmax(spart, good_SEBA)
    p
end

function _plot_SEBA(spart:: SEBAPartition{MultilayerGraph{T}, N}; kwargs...) where {T, N}
    h = []
    inds = 1:size(spart.vecs)[2]
    for i in inds
        push!(h,heatmap(reshape( spart.vecs[:,i],spart.partition.graph.N,spart.partition.graph.T),title="SEBA $(i), cut: $(round(spart.cuts[i], digits=3))",c=:tempo, xticks=get_ticks(spart.partition.graph.T), clims=(0,1), kwargs...))
    end
    h
end

function _plot_SEBA_trend(seba_part :: SEBAPartition{MultilayerGraph{T},N}, R) where {T <: Multiplex, N}
    cuts_multipleR = []
    for k in 1:R
        s_temp = SEBAPartition(seba_part.partition, k)
        push!(cuts_multipleR, s_temp.cuts)
    end

    cuts_full = zeros(2*R, R)
    for i in 1:R
        cuts_full[:,i] = [cuts_multipleR[i]; zeros(2*R-2*i)]
    end
    
    p_cuts = groupedbar(cuts_full', bar_position = :dodge, bar_width=0.7, label="", dpi=300,title=L"H(\mathcal{X}_k), \ k=1 \dots 2R", xlabel=L"R",titlefontsize=10)
    means = [sum(cuts_full[:,i])/(2*i) for i in 1:R]
    scatter!(p_cuts, means, label="")
    plot!(p_cuts, means, label="", c=:black)
    p_cuts
end

function _plot_SEBA_heatmapmax(spart :: SEBAPartition, good_SEBA :: Vector{Int64})
    N = spart.partition.graph.N
    T = spart.partition.graph.T
    max_vals, max_inds = findmax(spart.vecs[:,good_SEBA], dims=2)
    max_vals = reshape(max_vals, N,T)
    max_inds = reshape([x[2] for x in max_inds], N, T)
    p = nothing
    vals = []
    for i in 1:length(good_SEBA)
        vals_temp = zeros(N, T); vals_temp .= NaN
        vals_temp[findall(x->x==i, max_inds)] = max_vals[findall(x->x==i, max_inds)]
        if i == 1
            p = heatmap(vals_temp, c= cgrad([:white, color_vector(1)], [0.0,1.0]), dpi=300)
        else
            heatmap!(vals_temp, c= cgrad([:white, color_vector(i)], [0.0,1.0]), dpi=300)
        end

        push!(vals, vals_temp)
    end
   
    for i in 1:N, j in 1:T
        v = [vals[k][i,j] for k in 1:length(vals)]
        if sum(v .> 0.0) == 0
            max_inds[i,j] = 0
        end
    end

    p, vals, max_inds
end
