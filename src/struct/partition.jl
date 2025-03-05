"""
    SpectralPartition{MG <: MultilayerGraph, N <: Normalization}(mlgraph :: MultilayerGraph; compute_a = SpatTempMatching(), norm = IdentityNormalization())

Spectral partition instance for the temporal network defined by `mlgraph :: MultilayerGraph`. `mlgraph` is stored in `graph`, `norm :: Normalization` stores the Laplacian normalization function (default is `IdentityNormalization()`), `a` contains the diffusion constant, `L_temp` and `L_spat` are the temporal and spatial supra-Laplacians respectively and `evecs` and `evals` store the eigendata of the supra-Laplacian.

`SpectralPartition` computes the supra-Laplacian ``\\mathbf L^{(a)} = \\mathbf L^{\\rm spat} + a^2 \\mathbf L^{\\rm temp}``, computes the appropriate value of the diffusion constant ``a`` and computes the spectrum for both multiplex and nonmultiplex type networks.

Note that for nonmultiplex networks, `evecs` are lifted to ``\\mathbb R^{TN}``.
"""
mutable struct SpectralPartition{MG <: MultilayerGraph, N <: Normalization}
    graph :: MG
    norm :: N
    a :: Float64
    L_spat :: Union{Matrix{Float64}, SparseMatrixCSC}
    L_temp :: Union{Matrix{Float64}, SparseMatrixCSC}
    evecs :: Union{Matrix{Float64}, SparseMatrixCSC}
    evals :: Vector{Float64}
end

function Base.show(io::IO, ::MIME"text/plain", partition::SpectralPartition{MG,N}) where {MG,N}
    print("Spectral partition $(typeof(partition)) with $(typeof(partition.graph.connect)) connectivity and $(N) normalization. Diffusion constant a = $(partition.a).")
end

function _SpectralPartition(W_spat :: Matrix, W_temp :: Matrix, a, norm)
    infL = lap(W_spat) + a^2 .* lap(W_temp) |> norm
    evecs = eigvecs(infL)
    evals = eigvals(infL)
    return infL, evecs, evals
end

function _SpectralPartition(W_spat :: SparseMatrixCSC, W_temp :: SparseMatrixCSC, n_pairs,  a, norm)
    infL = lap(W_spat) + a^2 .* lap(W_temp) |> norm
    evals, evecs, sol = eigsolve(infL, rand(size(infL)[2]), n_pairs, :SR)
    @info sol
    return infL, hcat(evecs...), evals
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

function SpectralPartition(mlgraph :: MultilayerGraph{M}, a :: Float64, n_pairs :: Int64; norm = IdentityNormalization()) where M
    W_spat, W_temp = mlgraph.connect(mlgraph)
    infL, evecs, evals = _SpectralPartition(W_spat, W_temp, n_pairs,  a, norm)
    if M <: NonMultiplex
        evecs = _embed_compressed_evecs(mlgraph, evecs, Array(1:size(evecs)[2]), 0.0)
    end
    SpectralPartition(mlgraph, norm, a, lap(W_spat), lap(W_temp), evecs, evals)
end

function lift(partition :: SpectralPartition{MultilayerGraph{A}, B}) where {A <: NonMultiplex, B}
    return lift_matrix_nonmultiplex(partition.L_spat, partition.graph), lift_matrix_nonmultiplex(partition.L_temp, partition.graph)
end
