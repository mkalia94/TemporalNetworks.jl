# leidenalg via PyCall. Make sure before running julia, conda 
# env is set to graph, that has numpy, scipy, igraph and leidenalg
ENV["PYTHON"] = "python"
import .JuliaCommunity as juliac

function leiden_slice(pn :: SpectralPartition)
    ind_active = find_active(pn.graph)
    edges_W = [edgelist(pn.graph.W[i][ind_active[i], ind_active[i]]) for i in 1:pn.graph.T]
    network_compr = layer -> DataFrame(from = [x[1] for x in edges_W[layer]],
                                       to = [x[2] for x in edges_W[layer]],
                                       weight = [pn.graph.W[layer][ind_active[layer][x[1]], ind_active[layer][x[2]]] for x in edges_W[layer]]) 
    
    labelsSliceBySlice = zeros(pn.graph.N, pn.graph.T)

    for i in 1:pn.graph.T
        nodes_compr = DataFrame(id = Array(1:length(ind_active[i])),
                                label = ["$(x)" for x in Array(1:length(ind_active[i]))],
                                importance = ones(1:length(ind_active[i])))
        jc = juliac.JuliaCommunityInstance(network_compr(i), nodes = nodes_compr, node_label_field = "label", edge_weighted = true, task_series = "demo")
        juliac.set_method(jc, method = jc.methods.modularity)
        labels_leiden =juliac.leiden.find_partition(jc.igraph, juliac.leiden.ModularityVertexPartition, n_iterations=-1, max_comm_size=50)
        labelsSliceBySlice[ind_active[i],i] = labels_leiden.membership .+ 1
    end

    labelsSliceBySlice, stitch_slice(labelsSliceBySlice, pn.graph, pn.a)

end

function leiden_full(partition :: SpectralPartition{MultilayerGraph{A}, B}) where {A, B}
    if A <: NonMultiplex
        L_spat, L_temp = lift(partition)
        W_full = L_spat + partition.a^2 .* L_temp |> adj
        @show size(W_full)
    else
        W_full = partition.L_spat + partition.a^2 .* partition.L_temp |> adj
    end
    edges_W = edgelist(W_full)
    network = DataFrame(from = [x[1] for x in edges_W], to = [x[2] for x in edges_W], weight=[W_full[x...] for x in edges_W])
    nodes = DataFrame(id = Array(1:mlgraph.N*mlgraph.T), label = ["$(x)" for x in 1:mlgraph.N*mlgraph.T], importance = ones(mlgraph.N*mlgraph.T))
    jc = juliac.JuliaCommunityInstance(network, nodes=nodes, node_label_field = "label", edge_weighted=true, task_series = "demo")
    juliac.set_method(jc, method = jc.methods.modularity)
    #jc.quality = 0.9
    #jc.γ = 0.8
    juliac.discover_communities(jc)

    # labels_leiden =juliac.leiden.find_partition(jc.igraph, juliac.leiden.ModularityVertexPartition, n_iterations=-1)
    labels_full = jc.memberships
end

function stitch_slice(labels, mlgraph, a)
    N, T = size(labels)
    _, W_temp = mlgraph.connect(mlgraph)
    W_temp .*= a^2
    W_temp = lift_matrix_nonmultiplex(W_temp, mlgraph)
    labels_final = deepcopy(labels)
    layer1 = labels[:,1]
    K_max = maximum(labels)
    for i in 1:(mlgraph.T-1)
        K = maximum(layer1) |> Int64
        layer2 = labels[:, i+1]
        labels_final[:,i] = layer1

        cost = zeros(Int64(maximum(layer1)), Int64(maximum(layer2)))
        for j in 1:size(cost)[1], k in 1:size(cost)[2]
            vertices_j = findall(x->x==Float64(j), layer1)
            vertices_k = findall(x->x==Float64(k), layer2)
            cost[j,k] = sum(W_temp[vertices_j .+ (i-1)*mlgraph.N, vertices_k .+ (i)*mlgraph.N])
        end
        sol = bipartite_minimum_weight_edge_cover(-cost)
        label_new = deepcopy(layer2)

        if size(sol) == size(cost) # merging
            col_sums = sum(sol, dims=1)
            merge = findall(x->x>1, col_sums[1,:])
            matches = findall(x->x==1, col_sums[1,:])

            for ind in matches
                j = findfirst(x->x==1, sol[ind,:])
                vertices_j = findall(x->x==Float64(j), layer2)
                label_new[vertices_j] .= Float64(ind)
            end
            
            assigned = []
            for ind in merge
                j = findfirst(x->x==1, sol[ind,:])
                vertices_j = findall(x->x==Float64(j), layer2)
                if vertices_j ∉ assigned
                    label_new[vertices_j] .= Float64(ind)
                end
            end

        else
            sol = sol'
            row_sums = sum(sol, dims=2)
            split = findall(x->x>1, row_sums[:,1])
            matches = findall(x->x==1, row_sums[:,1])
            

            for ind in matches
                j = findfirst(x->x==1, sol[ind,:])
                vertices_j = findall(x->x==Float64(j), layer2)
                label_new[vertices_j] .= Float64(ind)
            end
            
            assigned = []

            for ind in split
                j = findall(x->x==1, sol[ind,:])
                _,j_min_ind = findmin(cost[ind,j])
                j_min = j[j_min_ind]
                ctr = 1
                for elem in j
                    vertices_elem = findall(x->x==Float64(elem), layer2)
                    if elem == j_min
                        label_new[vertices_elem] .= Float64(ind)
                    else
                        label_new[vertices_elem] .= K + ctr
                        @show K
                        ctr +=1
                    end
                end

                K_max = K_max + ctr
            end
        end


        layer1 = label_new
    end
    
    labels_final[:,end] = layer1
    labels_final
end







