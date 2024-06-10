abstract type IntraBlockConnectivity end
abstract type Block end


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


