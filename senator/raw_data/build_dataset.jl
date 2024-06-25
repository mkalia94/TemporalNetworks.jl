using Pkg; Pkg.activate("../../.")
using FileIO, JLD2
using LinearAlgebra, SparseArrays
using DataFrames, CSV

include("build_state.jl")

# If building for the first time with recently downloaded data,
# please filter it first to remove `bad' rollcalls, as per [Waugh 2009]. 
# This can be done via including the filter_df! function in gen_data_senator
# in build_senator.jl
include("build_senator.jl")

# MOST CRITICAL: Define congress range to construct network for. Note that
# it is generally quite time-consuming to process the whole dataset at once.
cong_range = 100:110

# Load CSV files
df = CSV.read("Sall_votes_filtered.csv", DataFrame)
n_cong = maximum(df.congress)
df2 = CSV.read("HSall_members.csv", DataFrame)

# Set up senator network data
voters, rolls, _, mat,  W_spat = gen_data_senator(df,cong_range)

# Set up state network data
_, W_spat_state, labels_statelist = gen_data_state(df, df2, voters, rolls, mat, cong_range)

W_spat_state = [Matrix{Float64}(x) for x in W_spat_state]

# Save state network
FileIO.save("../state.jld2", "W", W_spat_state, "labels", labels_statelist)


sens = []
for i in 1:length(cong_range)
    ind = findall(x->x>0, sum(W_spat[i], dims=1)[1,:])
    push!(sens, ind)
end

labels_party =  [df.party[findfirst(x->x==Int64(y), df.icpsr)] for y in voters[union(sens...)] ]

labels_state = [df2.state_abbrev[findfirst(x->x==Int64(y), df2.icpsr)] for y in voters[union(sens...)] ]

# Ensures nonmultiplex senator network is saved.
W_spat_trimmed = [Matrix{Float64}(W_spat[i][union(sens...), union(sens...)]) for i in 1:length(cong_range)]

# Save senator network data
FileIO.save("../senator.jld2","W", W_spat_trimmed, "labels_party", labels_party, "labels_state", labels_state)


