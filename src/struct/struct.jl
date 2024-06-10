abstract type Normalization end

"""
    (:: DegreeNormalization)(L :: Matrix)

Normalizes Laplacian with resepct to the standard normalization ``D^{-1/2}LD^{-1/2}`` where `D` is diagonal with `i`-th entry equal to `deg(i)`. 
"""
struct DegreeNormalization <: Normalization end

function (::DegreeNormalization)(L :: Matrix)
    sqrt.(deg(L))*L*sqrt.(deg(L))
end

"""
    `(:: IdentityNormalization)(L :: Matrix)`

Trivial Normalization of the Laplacian. Simply returns L. 
"""
struct IdentityNormalization <: Normalization end

function (::IdentityNormalization)(L :: Matrix)
    L
end

function deg(i :: Int64, L :: Matrix)
    return L[i,i]
end

function deg(L :: Matrix)
    return diagm(diag(L))
end

abstract type TemporalConnectivity end

"""
    (:: Multiplex)(mlgraph :: MultilayerGraph)

Connection type `:: TemporalConnectivity` for `MultilayerGraph`. Returns `Wˢᵖᵃᵗ, Wᵗᵉᵐᵖ` such that ``Wˢᵖᵃᵗ = ⨁ᵢWⁱ`` and ``Wᵗᵉᵐᵖ = W'⊗Iₙ``. Equivalently,

1. `W_spat = directsum(mlgraph.W)` where `W :: Vector{Matrix}`
2. `W_temp = kron(Wt, Matrix{Float64}(I, mlgraph.N, mlgraph.N))` where `Wt` is a `T x T` matrix describing how the layers are connected (symmetric supradiagonal).
"""
struct Multiplex <: TemporalConnectivity end

abstract type NonMultiplex <: TemporalConnectivity end

"""
    (:: NonMultiplexCompressed)(mlgraph :: MultilayerGraph)

Nonmultiplex connection type `:: TemporalConnectivity` for `MultilayerGraph`. Returns `W_spat, W_temp` such that `W_spat = directsum(mlgraph.W)` where `W :: Vector{Matrix}` and `W_temp[i,j] ≂̸ 0` iff `mlgraph.W[k][q,:] ≂̸ 0` and `mlgraph.W[l][r,:] ≂̸ 0` where ``(k,q) ∈ Rᴺ × Rᵀ`` and ``(l,r) ∈ Rᴺ × Rᵀ`` are the separated spacetime indices corresponding to ``i ∈ Rᴺᵀ`` and ``j ∈  Rᴺᵀ`` respectively.

"""
struct NonMultiplexCompressed <: NonMultiplex end

struct NonMultiplexDense <: NonMultiplex
    # Dense needs to be developed in time
end

abstract type DiffusionEstimator end


"""
    (:: SpatTempMatching)(::SpatTempMatching)(mlgraph :: MultilayerGraph{T}, norm :: Normalization, L_spat :: Matrix{Float64}, L_temp :: Matrix{Float64}) where T <: Multiplex

Computes the diffusion constant ``a`` by matching the second spatial eigenvalue (first non-trivial spatial eigenvalue) to the first temporal eigenvalue for the supra-Laplacian corresponding to a multiplex temporal network. The multiplex structure guarantees a unique smallest intersection point. See [FroylandKoltai2023] and [AtnipFroylandKoltai2024] for details on this heuristic. 

"""
struct SpatTempMatching <: DiffusionEstimator end

"""
    SpatTempMatchingNonMultiplex(evec_ind :: Int64, temp_ind :: Union{Int64, Nothng}, thresh :: Float64)

Creates `DiffusionEstimator` instance to compute the diffusion constant ``a`` for the supra-Laplacian arising from a non-multiplex network by matching the eigenvalue corresponding to `evec_ind` to the corresponding multiplex temporal eigenvalue with index `temp_ind`. The parameter `thresh` is used by `isspat(::MultilayerGraph{NonMultiplex}, ...)` to check for spatial-like behaviour.

    (:: SpatTempMatchingNonMultiplex)(::SpatTempMatching)(mlgraph :: MultilayerGraph{T}, norm :: Normalization, L_spat :: Matrix{Float64}, L_temp :: Matrix{Float64}) where T <: NonMultiplex

Computes the diffusion constant ``a`` using the matching heuristic for non-multiplex networks. See [FroylandKaliaKoltai2024] for details on this heuristic.

"""
struct SpatTempMatchingNonMultiplex <: DiffusionEstimator
    evec_ind :: Int64
    temp_ind :: Union{Int64, Nothing}
    thresh :: Float64
end

"""
    RayleighBalancing(ind :: Int64)

Creates `DiffusionEstimator` instance to compute the diffusion constant ``a`` for the supra-Laplacian arising from a non-multiplex network by matching the spatial and temproral Rayleigh quotient contributions with respect to the eigenvalue indexed by `ind`. Recommended for nonmultiplex networks. 

    (rb:: RayleighBalancing)(mlgraph :: MultilayerGraph ,norm :: Normalization, L_spat :: Matrix{Float64}, L_temp :: Matrix{Float64})


Computes the diffusion constant ``a`` using the Rayleigh balancing criterion. See [FroylandKaliaKoltai2024] for details on this heuristic.
"""
struct RayleighBalancing <: DiffusionEstimator
    ind :: Int64
end


