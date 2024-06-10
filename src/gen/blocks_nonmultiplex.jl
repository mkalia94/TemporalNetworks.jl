function (bll :: BlockGraphNonMultiplex)()
   # Make sure that the required cluster is the
   # first element of bll.clusters
    clusters, _, W_init = set_block_slice(bll.N, bll, 1)
    W_init[bll.clusters[1][3], bll.clusters[1][1]] .*= bll.interaction(1)
    W_init[bll.clusters[1][1], bll.clusters[1][3]] .*= bll.interaction(1)
    W_init[bll.clusters[1][2], :] .*= 0  
    W_init[:, bll.clusters[1][2]] .*= 0 
    W = []; push!(W, W_init)
    ctr = 0
    for i in 1:bll.T
        cluster_old = deepcopy(bll.clusters[1])
        bll.clusters[1][1] = [cluster_old[1][bll.evolve+1:end]; cluster_old[2][ctr + 1:ctr + bll.evolve]]
        bll.clusters[1][2] = [cluster_old[1][1:bll.evolve]; cluster_old[2][1:ctr]; cluster_old[2][ctr+bll.evolve+1:end-bll.evolve-ctr]; cluster_old[3][end-bll.evolve+1:end]; cluster_old[2][end-ctr+1:end]]
        bll.clusters[1][3] = [cluster_old[2][end-bll.evolve-ctr+1:end-ctr]; cluster_old[3][1:end-bll.evolve]]

        ctr += bll.evolve
        _,_,W_temp = set_block_slice(bll.N,bll,1)
        W_temp[bll.clusters[1][3], bll.clusters[1][1]] .*= bll.interaction(i)
        W_temp[bll.clusters[1][1], bll.clusters[1][3]] .*= bll.interaction(i)
        W_temp[bll.clusters[1][2], :] .*= 0  
        W_temp[:, bll.clusters[1][2]] .*= 0 
        push!(W,W_temp)
    end
    W
end


