"""
    MultilayerGraph{M <: TemporalConnectivity}(W :: Vector{Matrix{Float64}}; connect = Multiplex())

The temporal network instance. Stores sequence of adjacencies `W`, nodes per layer `N`, number of layers `T` and the connectivity of abstract type `TemporalConnectivity`. Default is `Multiplex()`, other options include `NonMultiplexCompressed()`. 
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
