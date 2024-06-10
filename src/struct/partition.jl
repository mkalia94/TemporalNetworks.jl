mutable struct SpectralPartition{MG <: MultilayerGraph, N <: Normalization}
    graph :: MG
    norm :: N
    a :: Float64
    L_spat :: Matrix{Float64}
    L_temp :: Matrix{Float64}
    evecs :: Matrix{Float64}
    evals :: Vector{Float64}
end

function Base.show(io::IO, ::MIME"text/plain", partition::SpectralPartition{MG,N}) where {MG,N}
    print("Spectral partition $(typeof(partition)) with $(typeof(partition.graph.connect)) connectivity and $(N) normalization. Diffusion constant a = $(partition.a).")
end

function _SpectralPartition(W_spat, W_temp, a, norm)
    infL = lap(W_spat) + a^2 .* lap(W_temp) |> norm
    evecs = eigvecs(infL)
    evals = eigvals(infL)
    return infL, evecs, evals
end

function SpectralPartition(mlgraph :: MultilayerGraph{M}; compute_a = SpatTempMatching(), norm = IdentityNormalization()) where M
    W_spat, W_temp = mlgraph.connect(mlgraph)
    a = compute_a(mlgraph, norm, lap(W_spat), lap(W_temp))
    infL, evecs, evals = _SpectralPartition(W_spat, W_temp, a, norm)
    if M <: NonMultiplex
        evecs = _embed_compressed_evecs(mlgraph, evecs, Array(1:size(evecs)[2]), 0.0)
    end
    SpectralPartition(mlgraph, norm, a, lap(W_spat), lap(W_temp), evecs, evals)
end

function SpectralPartition(mlgraph :: MultilayerGraph{M}, a :: Float64; norm = IdentityNormalization()) where M
    W_spat, W_temp = mlgraph.connect(mlgraph)
    infL, evecs, evals = _SpectralPartition(W_spat, W_temp, a, norm)
    if M <: NonMultiplex
        evecs = _embed_compressed_evecs(mlgraph, evecs, Array(1:size(evecs)[2]), 0.0)
    end
    SpectralPartition(mlgraph, norm, a, lap(W_spat), lap(W_temp), evecs, evals)
end

function lift(partition :: SpectralPartition{MultilayerGraph{A}, B}) where {A <: NonMultiplex, B}
    return lift_matrix_nonmultiplex(partition.L_spat, partition.graph), lift_matrix_nonmultiplex(partition.L_temp, partition.graph)
end
