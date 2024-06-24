function project_sens_to_state(partition :: SpectralPartition, evecs_temp :: Vector, labels_state :: Vector, order_states_ :: Vector)
    labelss = vcat([[x; x] for x in order_states_]...)  
    state_proj = [zeros(2*length(order_states_), partition.graph.T) for i in 1:length(evecs_temp)]
    sens = find_active(partition.graph)
    for i in 1:length(state_proj)
        for x in 1:length(order_states_), y in 1:size(state_proj[i])[2]
            ind = findall(z->(z ∈ sens[y] && labels_state[z] == order_states_[x]), 1:partition.graph.N)
            if length(ind)==1
                state_proj[i][(x-1)*2+1,y] = evecs_temp[i][ind[1],y]
            else
                state_proj[i][(x-1)*2+1:2*x,y] = evecs_temp[i][ind[1:2], y]
                # None of this accounts for states that may have more than two senators in a single congress (with incomplete terms)
            end
        end
    end
    # Filter
    for i in 1:length(state_proj)
        for (j,x) in enumerate(state_proj[i])
            if x == 0.0
                state_proj[i][j] = NaN
            end
        end
    end
    state_proj
end

