```@meta
CurrentModule = TemporalNetworks
```

# TemporalNetworks.jl 
`TemporalNetworks.jl` is a package that computes partitions of graphs/networks defined in spacetime. The networks are defined by a sequence `Vector{Matrix{Float64}}` of adjacnecy matrices. One can perform the following tasks:

- Compute spectral partitions of multiplex and nonmultiplex spacetime graphs using appropriate supra-Laplacians.
- Use the Sparse Eigenbasis Algorithm to disentagle multiple partition elements from eigenvectors of the supra-Laplacian.
- Plot networks, partitions and compute corresponding Cheeger ratios.
- Compute and plot partitions using Leiden, a modularity maximisation algorithm via the python package `leidenalg`.



