function gen_data_state(df_senator, df_state, voters_, rolls_, mat_, cong_range_)
    order_states = reverse(["HI","AK","WA", "OR", "CA", "WY", "UT", "NM", "NV", "MT", "ID", "CO", "AZ", "WV", "TN", "OK", "MD", "KY", "TX", "SC", "NC", "MS", "LA", "GA", "FL", "AR", "AL", "VA", "SD", "ND", "NE", "MO", "MN", "KS", "IA", "WI", "OH", "MI", "IN", "IL", "PA", "NY", "NJ", "DE", "VT", "RI", "NH", "MA", "ME", "CT", ""])

    dff = innerjoin(df_senator,df_state, on=[:icpsr, :congress, :chamber])
    states = unique(dff.state_abbrev)
    if length(states) == 50
        states = order_states[2:end]
    end

    unique_states = []
    for voter in voters_
        push!(unique_states,unique(dff.state_abbrev[findall(x->x==voter, dff.icpsr)]))
    end

    ind_states = zeros(length(voters_))

    for i in 1:length(unique_states)
        try
            ind = findfirst(x->x==unique_states[i][1], states)
            ind_states[i] = ind
        catch e
            @show e
        end
    end

    mat_state = zeros(length(states), sum(rolls_.counts))

    for i in 1:length(ind_states)
        if ind_states[i] != 0.0
            mat_state[Int64(ind_states[i]),:] += mat_[i,:]
        end
    end

    mat_split = []

    for i in cong_range_ #  1:maximum(df_senator.congress)
        ind = sum(rolls_.counts[1:(i-1)]) + 1: sum(rolls_.counts[1:i])
        push!(mat_split, mat_state[:,ind])
    end

    M = length(cong_range_)
    N = size(mat_split[1])[1]

    # Correlation based (perhaps useful for signed networks)
    # W_spat_state = [cor(mat_split[i]') for i in 1:M]
    # for i in 1:M
    #     W_spat_state[i][diagind(W_spat_state[i])] .= 0.0
    # end

    mat_split, W_spat_state = construct_mat(df = dff, mat = mat_state, object = states, rolls_ = rolls_, cong_range_ = cong_range_)
    mat_split, W_spat_state, states
end

function construct_mat(; df, mat, object, rolls_, cong_range_)
    mat_split = []
    for i in cong_range_ #  1:maximum(df.congress)
        ind = sum(rolls_.counts[1:(i-1)]) + 1: sum(rolls_.counts[1:i])
        push!(mat_split, mat[:,ind])
    end
    
    W_spat = []
    for Wraw in mat_split
        Wtemp = sparse((abs.(Wraw)*abs.(Wraw)'))
        W = sparse(zeros(length(object), length(object)))
        x,y,_ = findnz(Wtemp)
        for i in 1:length(x)
            W[x[i],y[i]] = sum((Wraw[x[i],:] .* Wraw[y[i],:]) .> 0)/Wtemp[x[i],y[i]]
        end  
        W[diagind(W)] .= 0
        push!(W_spat, dropzeros(W))
    end
    mat_split, W_spat
end

