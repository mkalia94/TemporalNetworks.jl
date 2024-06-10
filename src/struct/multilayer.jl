"""
    MultilayerGraph{M <: TemporalConnectivity}(W :: Vector{Matrix{Float64}}, N :: Union{Int64, Float64}, T :: Int64)

Creates `DiffusionEstimator` instance to compute the diffusion constant ``a`` for the supra-Laplacian arising from a non-multiplex network by matching the spatial and temproral Rayleigh quotient contributions with respect to the eigenvalue indexed by `ind`. Recommended for nonmultiplex networks. 

    (rb:: RayleighBalancing)(mlgraph :: MultilayerGraph ,norm :: Normalization, L_spat :: Matrix{Float64}, L_temp :: Matrix{Float64})
`

Computes the diffusion constant ``a`` using the Rayleigh balancing criterion. See [FroylandKaliaKoltai2024] for details on this heuristic.
"""

mutable struct MultilayerGraph{M <: TemporalConnectivity}
    W :: Vector{Matrix{Float64}}
    N :: Union{Int64, Float64}
    T :: Int64
    connect :: M
end

function Base.show(io::IO, ::MIME"text/plain", mlgraph::MultilayerGraph{M}) where M
    print("Multilayer graph with $(M) connectivity. $(mlgraph.T) time slices with $(mlgraph.N) nodes each.")
end
    

function MultilayerGraph(W; connect = Multiplex())
    N = [size(x)[1] for x in W]
    if length(unique(N)) == 1
        return MultilayerGraph(W, only(unique(N)), length(W), connect)
    else
        return MultilayerGraph(W, N, length(W), connect)
    end
end

function (ff::Multiplex)(mlgraph :: MultilayerGraph)
    Wt = zeros(mlgraph.T,mlgraph.T)
    for i in 1:mlgraph.T, j in 1:mlgraph.T
        if i == j+1 || i == j-1
            Wt[i,j] = 1.0
        end
    end
    directsum(mlgraph.W), kron(Wt, Matrix{Float64}(I, mlgraph.N, mlgraph.N))
end


function _Layer2Layer(T)
    Wt = zeros(T,T)
    for i in 1:T, j in 1:T
        if i == j+1 || i == j-1
            Wt[i,j] = 1.0
        end
    end
    Wt
end
