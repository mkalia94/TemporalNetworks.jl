function (::NonMultiplexCompressed)(mlgraph)
    ind_active = find_active(mlgraph)
    W_spat = [mlgraph.W[i][ind_active[i], ind_active[i]] for i in 1:mlgraph.T]
    W_spat_full = directsum(typeof(mlgraph.W)(W_spat))

    # Wtemp
    ctr = 0
    ind_active_full = deepcopy(ind_active)
    for (i, vec) in enumerate(ind_active)
        ind_active_full[i] = vec .+ ctr
        ctr += mlgraph.N
    end

    ind_active_full = vcat(ind_active_full...)
    W_temp_full = Multiplex()(mlgraph)[2][ind_active_full, ind_active_full]
    
    return W_spat_full, W_temp_full
end

function find_active(mlgraph:: MultilayerGraph)
    ind_active = []
    for i in 1:mlgraph.T
        ind = findall(x->x>0, sum(mlgraph.W[i], dims=2)[:,1])
        push!(ind_active, ind)
    end
    ind_active
end

function embed_compressed_evecs(partition :: SpectralPartition{MultilayerGraph{N}, M}, inds) where {N <: NonMultiplex, M <: Normalization}
    _embed_compressed_evecs(partition.graph, partition.evecs, inds, NaN)
end

function embed_compressed_evecs(spart :: SEBAPartition{MultilayerGraph{N}, M}) where {N <: NonMultiplex, M <: Normalization}
    _embed_compressed_evecs(spart.partition.graph, spart.vecs, 1:length(spart.inds), 0.0)
end

function _embed_compressed_evecs(graph :: MultilayerGraph{M}, vecs :: Matrix{Float64}, inds :: Vector{Int64}, factor :: Float64) where M
    ind_active = find_active(graph)
    evecs_embed = zeros(graph.N*graph.T, length(inds))
    for (k, ind) in enumerate(inds)
        evec_full = zeros(graph.N,graph.T) .* factor
        ctr = 1
        for i in 1:graph.T
            evec_full[ind_active[i],i] = vecs[ctr:ctr+length(ind_active[i])-1,ind]
            ctr = ctr + length(ind_active[i])
        end
        evecs_embed[:,k] = reshape(evec_full, graph.N*graph.T)
    end

    evecs_embed
end

function lift_matrix_nonmultiplex(mat :: Union{Matrix, SparseMatrixCSC}, mlgraph :: MultilayerGraph)
    inds_active = find_active(mlgraph)
    inds_active = vcat([x .+ (i-1)*mlgraph.N for (i,x) in enumerate(inds_active)]...)
    
    if typeof(mat) <: Matrix
        lift_mat = zeros(mlgraph.N*mlgraph.T, mlgraph.N*mlgraph.T)
    else
        lift_mat = spzeros(mlgraph.N*mlgraph.T, mlgraph.N*mlgraph.T)
    end
    lift_mat[inds_active, inds_active] = mat
    lift_mat
end

function _isspat_nonmultiplex(mlgraph :: MultilayerGraph{A}, vec; indx = nothing) where A <: NonMultiplex
    W_temp = _Layer2Layer(mlgraph.T)
    evec_temp = eigvecs(lap(W_temp))
    evec_temp_perfect = [kron(evec_temp[:,i], ones(mlgraph.N)) for i in 1:mlgraph.T]
    evec_temp_perfect = hcat([x ./ sqrt(sum(abs2, x)) for x in evec_temp_perfect]...)
    if isnothing(indx)
        _, ii = findmax([abs(sum(vec .* evec_temp_perfect[:,i])) for i in 1:mlgraph.T])
        return sum(vec .* evec_temp_perfect[:,ii])
    else
        return sum(vec .* evec_temp_perfect[:, indx])
    end

end

function isspat(mlgraph :: MultilayerGraph{A}, norm :: Normalization,  L_spat, L_temp, a ::  Float64, ind :: Int64, indx :: Union{Nothing, Int64}) where A <: NonMultiplex
    evecs = eigvecs(L_spat + a^2 .* L_temp)
    evecs = _embed_compressed_evecs(mlgraph, evecs, Array(1:size(evecs)[2]), 0.0)
    _isspat_nonmultiplex(mlgraph, evecs[:,ind], indx = indx)
end
