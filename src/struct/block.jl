abstract type IntraBlockConnectivity end
abstract type Block end

"""
    BlockGraph(N :: Int64, T :: Int64, list :: Vector, η :: Float64, clusters :: Union{Nothing, Vector}, degrees :: Union{Nothing, Vector})

Multiplex type block graph instance. Constructs a temporal network of `N` vertices per time slice with `T` time steps.

The transition types are encoded in `list` which is a `Vector{Int64}`. Numbers indicate number of disjoint clusters to be constructed at certain time points. `clusters` and `degrees` hold information about which vertices belong to which clusters and the cluster degrees respectively (the clusters are regular subgraphs). 

The instance used as a function and passed with no arguments returns a sequence of adjacency matrices of type `Vector{Float64}`. 
"""
mutable struct BlockGraph <: Block
    N :: Int64 
    T :: Int64
    list :: Vector
    η :: Float64 # fraction of nodes within a community at t=0
    clusters :: Union{Nothing, Vector}
    degrees :: Union{Nothing, Vector}
end

struct ScalingConnectivity <: IntraBlockConnectivity
    factor :: Float64
end

function (sc :: ScalingConnectivity)(x :: Int64)
    x/sc.factor
end

"""
    BlockGraphNonMultiplex(N :: Int64, T :: Int64, list :: Vector, η :: Float64, clusters :: Union{Nothing, Vector}, degrees :: Union{Nothing, Vector}, interaction :: IntraBlockConnectivity, evolve :: Int64)

Non-multiplex type staircase graph instance. Constructs a temporal network of `N`spatial vertices `\bigcup_{t=1}^T\{x: (t,x)\in\mathcal{V}\}`` with `T` time steps. 

The transition types are encoded in `list` which is a `Vector{Int64}`. Numbers indicate number of disjoint clusters to be constructed at certain time points. `clusters` and `degrees` hold information about which vertices belong to which clusters and the cluster degrees respectively (the clusters are regular subgraphs). 

`interaction` is a `IntraBlockConnectivity` instance that acts as a function, and is applied to all intercluster edge weights. The default is `ScalingConnectivity(T)` which is the function ``x \mapsto x/T``.  

`evolve` indicates how many vertices per block are to be considered absent in the next layer, and switched with equally many new vertices from the absent ones. Defaults to 1.


    BlockGraphNonMultiplex(N, T, list, η, clusters, degrees)

`BlockGraphNonMultiplex` instance with `interactivity = ScalingConnectivity(2)` and `evolve=1`.

    BlockGraphNonMultiplex(N, T, list, η, clusters, degrees)

`BlockGraphNonMultiplex` instance with `interactivity = ScalingConnectivity(2)`

The instance used as a function and passed with no arguments returns a sequence of adjacency matrices of type `Vector{Float64}`. 
"""
struct BlockGraphNonMultiplex <: Block
    N :: Int64 
    T :: Int64
    list :: Vector
    η :: Float64 # fraction of nodes within a community at t=0
    clusters :: Union{Nothing, Vector}
    degrees :: Union{Nothing, Vector}
    interaction :: IntraBlockConnectivity
    evolve :: Int64
end
    
function BlockGraphNonMultiplex(N, T, list, η, clusters, degrees)
    BlockGraphNonMultiplex(N, T, list, η, clusters, degrees, ScalingConnectivity(T), 1)
end

function BlockGraphNonMultiplex(N, T, list, η, clusters, degrees, evolve)
    BlockGraphNonMultiplex(N, T, list, η, clusters, degrees, ScalingConnectivity(T), evolve)
end


