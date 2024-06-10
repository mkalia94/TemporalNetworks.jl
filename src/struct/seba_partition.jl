mutable struct SEBAPartition{MG <: MultilayerGraph, N <: Normalization}
    partition :: SpectralPartition{MG, N}
    inds :: Vector{Int64}
    vecs :: Matrix{Float64}
    cuts :: Vector{Float64}
end

function Base.show(io::IO, ::MIME"text/plain", s_partition::SEBAPartition{MG,N}) where {MG,N}
    print("SEBA partition $(typeof(s_partition)) with $(size(s_partition.vecs)[2]) vectors using $(["$(x), " for x in s_partition.inds[1:end-1]]...)$(s_partition.inds[end]) spatial eigenvectors.")
end

function Base.length(s_partition :: SEBAPartition)
    size(s_partition.vecs)[2]
end

function SEBAPartition(partition :: SpectralPartition{MultilayerGraph{N}, M}, inds :: Vector{Int64}) where {N <: TemporalConnectivity, M <: Normalization}
    vecs, rot = SEBA(double_basis_SEBA(partition.graph, inds, partition.evecs)) 
    cuts = cuts_SEBA(partition, vecs) 
    SEBAPartition(partition, inds, vecs, cuts)
end

function SEBAPartition(partition :: SpectralPartition{MultilayerGraph{N}, M}, max_ind :: Int64) where {N <: Multiplex, M <: Normalization}

    evecs_for_SEBA = []
    for i in 2:15
        if isspat(partition.evecs[:,i], partition.graph.N, partition.graph.T, 0.01) # This doesn't work for non-multiplex type networks
            push!(evecs_for_SEBA,i)
        end
        length(evecs_for_SEBA) >= max_ind ? break : 0
    end
    SEBAPartition(partition, convert(Vector{Int64},evecs_for_SEBA))
end

function SEBAPartition(partition :: SpectralPartition{MultilayerGraph{N},M}) where {N <: Multiplex, M <: Normalization}
    SEBAPartition(partition, 3)
end

function SEBAPartition(partition :: SpectralPartition{MultilayerGraph{N},M}) where {N <: NonMultiplex, M <: Normalization}
    SEBAPartition(partition, [1,2,3])
end

