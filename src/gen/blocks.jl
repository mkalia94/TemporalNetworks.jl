function (bl :: BlockGraph)()  
    
    W = []
    clusters_full = []

    cluster1,_,W1 = set_block_slice(bl.N, bl, 1)
    cluster2 = nothing
    W2 = nothing
    for i in 1:length(bl.list)-1
        cluster2,_,W2 = set_block_slice(bl.N, bl, i+1)
        if bl.list[i+1] == 0
            W_merge = reverse(blend_slices(W2, cluster2, W1, cluster1, bl.T))
        else
            W_merge = blend_slices(W1, cluster1, W2, cluster2, bl.T)
        end
        push!(W, W_merge...)
        push!(clusters_full, cluster1)
        if i == (length(bl.list)-1)
            push!(clusters_full, cluster2)
        end
        cluster1 = deepcopy(cluster2)
        W1 = deepcopy(W2)
    end
    bl.clusters = clusters_full
    Vector{Matrix{Float64}}(W)
end

function blend_slices(first_, clusters_first, second_, clusters_second, T)
    
    half_way = T/2

    W = []

    all_edges = union!(edgelist(first_),edgelist(second_))
    edges_delete = [z for z ∈ edgelist(first_) if z ∉ edgelist(second_) ]
    edges_add = [z for z ∈ edgelist(second_) if z ∉ edgelist(first_) ]
    n_actions = length(edges_delete) + length(edges_add)
    actions_per = ceil(n_actions/T)   |> Int

    actions = []

    W_temp = deepcopy(first_)
    W_old = deepcopy(W_temp)
    push!(W,deepcopy(W_temp))
    layer = 1

    while W_temp != second_# length(actions) < n_actions
        sample_action = shuffle(all_edges)[1:actions_per]
        for edge in sample_action
            if edge ∉ actions && edge ∈ edges_add
                W_temp[edge...] = 1.0
                W_temp[reverse(edge)...] = 1.0
                push!(actions,edge)
            elseif edge ∉ actions && edge ∈ edges_delete
                W_temp[edge...] = 0.0
                W_temp[reverse(edge)...] = 0.0
                push!(actions,edge)
            end
        end
        if W_temp != W_old && ((layer < half_way && minimum(comm_quality(W_temp,clusters_first)) >=0.6) || (layer >= half_way && minimum(comm_quality(W_temp,clusters_second))>=0.6)) 
            push!(W,deepcopy(W_temp))
            W_old = deepcopy(W_temp)
            layer = layer + 1
        end
    end

    W
end


function set_block_slice(N, bl :: A, ind_block) where A <: Block
    if isnothing(bl.clusters)
        return set_block_slice_fromscratch(N, bl, ind_block)
    else
        return set_block_slice_apriori(N, bl, ind_block)
    end
end

function set_block_slice_apriori(N, bl :: A, ind_block) where A <: Block
    W = zeros(N,N)
    clusters = bl.clusters[ind_block]
    degs = bl.degrees[ind_block]
    ctr = 1
    for cluster in clusters
        W[cluster, cluster] = partiallyconnected(length(cluster), degs[ctr])
        ctr = ctr + 1
    end

    if bl.list[ind_block] > 0
        edge_per = floor((1-bl.η)*minimum(length.(clusters)))
        edge_per = minimum([floor((1-bl.η)*minimum(degs)),edge_per]) |> Int
        nodes_extraedges = [shuffle(clusters[i])[1:edge_per] for i in 1:(bl.list[ind_block]+1)]

        extra_edges = make_edges(nodes_extraedges)

        for edge in extra_edges
            W[edge...] = 1.0
            W[reverse(edge)...] = 1.0
        end
    else
        nodes_extraedges = nothing
    end

    return clusters, nodes_extraedges, W

end

function set_block_slice_fromscratch(N,bl :: A, ind_block) where A <: Block

    n = bl.list[ind_block]
    η = bl.η

    # 1. Get the right sizes
    size_ = floor(N/(n+1.5)) |> Int
    max_degree = minimum([size_, N - n*size_])
    max_degree = 4
    if isodd(N - n*size_) && isodd(max_degree)
        max_degree = max_degree + 1
    end
    if n == 0
        max_degree = floor(N/2)+1 |> Int
        max_degree = 4
    end

    println("------- $(n) clusters of size $(size_) each with remaining $(N-n*size_) nodes of degree $(max_degree) -------- ")

    # 2. Project clusters onto labelled nodes
    nodes = Array(1:N)
    clusters = isnothing(bl.clusters) ? [] : bl.clusters[ind_block]
    if n > 0
        size_rest = N - n*size_
        rest_ind = shuffle(nodes[1:N - size_rest - 1])[1]
        rest_ind = Array(rest_ind:(rest_ind + size_rest - 1))
        push!(clusters, nodes[rest_ind])
        println("Cluster 1: $(clusters[1])")
        comm_ind = mod1(rest_ind[end] + 1, N)
        for i in 1:n
            ind = mod1.(Array(comm_ind: comm_ind + size_-1),N)
            push!(clusters, nodes[ind])
            println("Cluster $(i+1): $(clusters[i+1])")
            comm_ind = comm_ind + size_
        end
        W = zeros(N,N)
        ctr = 1
        for cluster in clusters
            if ctr == 1
                W[cluster, cluster] = partiallyconnected(size_rest, max_degree)
            else
                W[cluster, cluster] = fullyconnected(size_)
            end
            ctr = ctr + 1
        end

    else
        push!(clusters, Array(1:N))
        W = partiallyconnected(N, max_degree)
    end

    # 3. Stitch the disjoint clusters (connect communities to the rest)
    nodes_extraedges = []
    if n > 0
        edge_per = floor((1-η)*size_)
        edge_per = minimum([floor((1-η)*max_degree),edge_per]) |> Int
        nodes_extraedges = [shuffle(clusters[i])[1:edge_per] for i in 1:(n+1)]

        extra_edges = make_edges(nodes_extraedges)

        for edge in extra_edges
            W[edge...] = 1.0
            W[reverse(edge)...] = 1.0
        end
    end

    return clusters, nodes_extraedges, W
end


function make_edges(clusters)
    
    # Find rest cluster
    rest = clusters[1]

    edges = []
    for i in 1:length(clusters)
        if i == 1
            continue
        end
        for j in 1:length(clusters[i])
            # Choose new vertex from another cluster
            for v in rest
                edge = [clusters[i][j], v]
                if edge ∉ edges && reverse(edge) ∉ edges
                    push!(edges, edge)
                end
            end
        end
    end
    edges
end

function edgelist(W)
    edges = []
    for i in 1: size(W)[1], j in 1:size(W)[2]
        if W[i,j] > 0
            push!(edges,[i,j])
        end
    end
    edges
end

function comm_quality(W :: Matrix, V :: Vector)
    [comm_quality(W, cluster) for cluster in V]
end

function comm_quality(W:: Matrix, V::Vector{Int64})
    edges_within = sum(W[edge...] for edge in collect(product(V,V)))
    edges_total = sum(W[edge...] for edge in collect(product(V, Array(1:size(W)[1]))))
    edges_within/edges_total
end

function partiallyconnected(N,d)
    W = zeros(N,N)
    if iseven(d)
        W = pc(N,d)
    elseif isodd(d) && isodd(N)
        error("Handshaking lemma: N and d should not both be odd")
    else
        W = pc(N,d-1)
        for i in 1:N
            v1 = i
            v2 = mod(v1-1 + Int(N/2),N)+1
            W[v1,v2] = 1.0
        end
    end
    W
end

function pc(N::Int64,d::Int64)
    if !iseven(d)
        error("function pc should be supplied with even d")
    end
    W = zeros(N,N)
    arr = [Array(-Int(d/2):-1); Array(1:Int(d/2))]
    for i in 1:N
        v1 = i
        for k in arr
            v2 = mod(v1-1+k,N)+1
            W[v1, v2] = 1.0
        end
    end
    W
end

function fullyconnected(N)
    return ones(N,N) .- diagm(N,N,ones(N))
end


