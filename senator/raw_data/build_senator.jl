function gen_data_senator(df_senator, cong_range)
    
    voters = unique(df_senator.icpsr)
    
    # Coment out in case the data is not filtered.
    # Here we remove `bad' roll calls.
    # filter_df!(df_senator)

    transform!(df_senator, :icpsr => ByRow(x->findfirst(y->y==x, voters)) => :voter)
    gdf2 = groupby(df_senator, [:congress])
    rolls = combine(gdf2, :rollnumber => maximum => :counts)
    gdf2 = nothing
    transform!(df_senator, [:congress, :rollnumber] => ByRow((x,y)->sum(rolls.counts[1:(x-1)])+y) => :ind)
    gdf = groupby(df_senator, [:voter])
    mat = sparse(zeros(length(voters), sum(rolls.counts)))
    for g in gdf
        mat[g.voter[1], g.ind] = g.cast_code
    end
    mat_split = []


    for i in cong_range #  1:maximum(df.congress)
        ind = sum(rolls.counts[1:(i-1)]) + 1: sum(rolls.counts[1:i])
        push!(mat_split, mat[:,ind])
    end
    W_spat = []; 
    for Wraw in mat_split
        Wtemp = sparse((abs.(Wraw)*abs.(Wraw)'))
        W = sparse(zeros(length(voters), length(voters)))
        x,y,_ = findnz(Wtemp)
        for i in 1:length(x)
            W[x[i],y[i]] = sum((Wraw[x[i],:] .* Wraw[y[i],:]) .> 0)/Wtemp[x[i],y[i]]
        end  
        W[diagind(W)] .= 0
        push!(W_spat, dropzeros(W))
    end
    
    voters, rolls, mat_split, mat,  W_spat


end

function filter_df!(df)
    for i in 1:size(df)[1]
        if df.cast_code[i] == 6.0
            df.cast_code[i] = -1.0
        elseif df.cast_code[i] != 1.0
            df.cast_code[i] = 0.0
        end
    end
    
    # Filter bad rollnumbers
    df_groups = groupby(df,[:congress, :rollnumber])
    df_stats = combine(df_groups, :cast_code => x -> abs(sum(x)/sum(abs,x)))
    bad_rolls = []
    for i in 1:size(df_stats)[1]
        if df_stats.cast_code_function[i] >= 0.94
            push!(bad_rolls, [df_stats.congress[i], df_stats.rollnumber[i]])
        end
    end
    
    
    filter!([:congress,:rollnumber]=> (x,y)->[x,y] ∉ bad_rolls, df)
end

