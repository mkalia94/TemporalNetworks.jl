using TemporalNetworks, LinearAlgebra
using Test

@testset "Multiplex case" begin
    # Set up toy model
    W1 = TemporalNetworks.fullyconnected(5)
    W5 = deepcopy(W1)
    W2 = deepcopy(W1)
    W2[[3,4],[4,3]] .= 0
    W2[[2,4],[4,2]] .= 0
    W2[[1,4],[4,1]] .= 0
    W2[[2,5],[5,2]] .= 0
    W2[[1,5],[5,1]] .= 0
    W3 = deepcopy(W2)
    W3[[3,5],[5,3]] .= 0
    W4 = deepcopy(W3)
    W4[[1,4],[4,1]] .= 1
    Wtoy = [W1, W2, W3, W4, W1]

    # Build temporal networks, spectral partitions and SEBA partitions
    mlgraph = MultilayerGraph(Wtoy)
    partition = SpectralPartition(mlgraph)
    spart = SEBAPartition(partition,1)
    Wt = zeros(5,5)
    ind = diag([(i-1)*5+j for i in 2:5, j in 1:4])
    Wt[ind] .= 1.0
    Wt = Wt + Wt'

    # Tests
    
    # Test connectivity
    @test mlgraph.connect(mlgraph) == (TemporalNetworks.directsum(Wtoy), kron(Wt, Matrix{Float64}(I,5,5)))

    # Test a-value for the unnormalized case
    @test partition.a ≈ 2 atol=0.05

    # Test the support of SEBA vectors for the correct partitions
    ind1 = sort(findall(x->x>0.05, spart.vecs[:,1]))
    ind2 = sort(findall(x->x>0.05, spart.vecs[:,2]))
    ind_gt1 = [9,10,14,15,19,20]
    ind_gt2 = [x for x in 6:20 if x∉ ind_gt1]
    @test ind1 == ind_gt1 || ind1 == ind_gt2
    @test ind2 == ind_gt1 || ind2 == ind_gt2 

    # Test cut values for the unnormalized case
    @test sort(spart.cuts) ≈ [2.29, 3] atol=0.01
    
    # Normalized case
    partition_norm = SpectralPartition(mlgraph, norm = DegreeNormalization())
    @test partition_norm.a ≈ 1.97 atol=0.01
    spart_norm = SEBAPartition(partition_norm,1)
    ind1 = sort(findall(x->x>0.05, spart_norm.vecs[:,1]))
    ind2 = sort(findall(x->x>0.05, spart_norm.vecs[:,2]))
    ind_gt1 = [9,10,14,15,19,20]
    ind_gt2 = [x for x in 6:20 if x∉ ind_gt1]
    @test ind1 == ind_gt1 || ind1 == ind_gt2
    @test ind2 == ind_gt1 || ind2 == ind_gt2
    @test sort(spart_norm.cuts) ≈ [0.281, 0.321] atol=0.01
end


@testset "Nonmultiplex case" begin
    
    η = 0.8
    list = [2,1]
    clusters = [[Array(1:2), Array(3:6), Array(7:8)],
                    [Array(1:10), Array(11:20)]] # you can ignore this line
    degrees = [[1,3,1],
            [6,6]] # Ignore this line
    evolve = 1
    block = BlockGraphNonMultiplex(8, 2, list, η, clusters, degrees, evolve)
    W2 = block() |> Vector{Matrix{Float64}}
    # Introduce some inter cluster connections
    W2[1][2,7] = 0.1; W2[1][7,2] = 0.1
    W2[2][3,6] = 0.1; W2[2][6,3] = 0.1
    W2[3][4,5] = 0.1; W2[3][5,4] = 0.1
    
    mlgraph_nm = MultilayerGraph(W2, connect = NonMultiplexCompressed())
    partition_nm = SpectralPartition(mlgraph_nm, compute_a = RayleighBalancing(2))
    spart_nm = SEBAPartition(partition_nm,[2,])
    ind1 = sort(findall(x->x>0, spart_nm.vecs[:,1]))
    ind2 = sort(findall(x->x>0, spart_nm.vecs[:,2]))
    ind_gt1 = [7,8,14,15,21,22]
    ind_gt2 = [1,2,10,11,19,20]
    @test ind1 == ind_gt1 || ind1 == ind_gt2
    @test ind2 == ind_gt1 || ind2 == ind_gt2
    @test partition_nm.a ≈ 0.472 atol=0.01
    
end

