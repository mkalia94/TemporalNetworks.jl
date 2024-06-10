abstract type Toy end

struct ToyMultiplex <: Toy
    N :: Int64
    T :: Int64
end

struct ToyStaircase <: Toy

end

function (toy::ToyMultiplex)()
    list = [0,1,2]
    degrees = nothing
    η = 0.8
    clusters = nothing
    block = BlockGraph(toy.N, toy.T, list, η, clusters, degrees)
    block()
end


