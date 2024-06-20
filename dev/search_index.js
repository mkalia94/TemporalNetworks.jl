var documenterSearchIndex = {"docs":
[{"location":"","page":"TemporalNetworks","title":"TemporalNetworks","text":"CurrentModule = TemporalNetworks","category":"page"},{"location":"#TemporalNetworks","page":"TemporalNetworks","title":"TemporalNetworks","text":"","category":"section"},{"location":"","page":"TemporalNetworks","title":"TemporalNetworks","text":"","category":"page"},{"location":"","page":"TemporalNetworks","title":"TemporalNetworks","text":"Modules = [TemporalNetworks]","category":"page"},{"location":"#TemporalNetworks.DegreeNormalization","page":"TemporalNetworks","title":"TemporalNetworks.DegreeNormalization","text":"(:: DegreeNormalization)(L :: Matrix)\n\nNormalizes Laplacian with resepct to the standard normalization D^-12LD^-12 where D is diagonal with i-th entry equal to deg(i). \n\n\n\n\n\n","category":"type"},{"location":"#TemporalNetworks.IdentityNormalization","page":"TemporalNetworks","title":"TemporalNetworks.IdentityNormalization","text":"`(:: IdentityNormalization)(L :: Matrix)`\n\nTrivial Normalization of the Laplacian. Simply returns L. \n\n\n\n\n\n","category":"type"},{"location":"#TemporalNetworks.Multiplex","page":"TemporalNetworks","title":"TemporalNetworks.Multiplex","text":"(:: Multiplex)(mlgraph :: MultilayerGraph)\n\nConnection type :: TemporalConnectivity for MultilayerGraph. Returns Wˢᵖᵃᵗ, Wᵗᵉᵐᵖ such that Wˢᵖᵃᵗ = ᵢWⁱ and Wᵗᵉᵐᵖ = WIₙ. Equivalently,\n\nW_spat = directsum(mlgraph.W) where W :: Vector{Matrix}\nW_temp = kron(Wt, Matrix{Float64}(I, mlgraph.N, mlgraph.N)) where Wt is a T x T matrix describing how the layers are connected (symmetric supradiagonal).\n\n\n\n\n\n","category":"type"},{"location":"#TemporalNetworks.NonMultiplexCompressed","page":"TemporalNetworks","title":"TemporalNetworks.NonMultiplexCompressed","text":"(:: NonMultiplexCompressed)(mlgraph :: MultilayerGraph)\n\nNonmultiplex connection type :: TemporalConnectivity for MultilayerGraph. Returns W_spat, W_temp such that W_spat = directsum(mlgraph.W) where W :: Vector{Matrix} and W_temp[i,j] ≂̸ 0 iff mlgraph.W[k][q,:] ≂̸ 0 and mlgraph.W[l][r,:] ≂̸ 0 where (kq)  Rᴺ  Rᵀ and (lr)  Rᴺ  Rᵀ are the separated spacetime indices corresponding to i  Rᴺᵀ and j   Rᴺᵀ respectively.\n\n\n\n\n\n","category":"type"},{"location":"#TemporalNetworks.RayleighBalancing","page":"TemporalNetworks","title":"TemporalNetworks.RayleighBalancing","text":"RayleighBalancing(ind :: Int64)\n\nCreates DiffusionEstimator instance to compute the diffusion constant a for the supra-Laplacian arising from a non-multiplex network by matching the spatial and temproral Rayleigh quotient contributions with respect to the eigenvalue indexed by ind. Recommended for nonmultiplex networks. \n\n(rb:: RayleighBalancing)(mlgraph :: MultilayerGraph ,norm :: Normalization, L_spat :: Matrix{Float64}, L_temp :: Matrix{Float64})\n\nComputes the diffusion constant a using the Rayleigh balancing criterion. See [FroylandKaliaKoltai2024] for details on this heuristic.\n\n\n\n\n\n","category":"type"},{"location":"#TemporalNetworks.SpatTempMatching","page":"TemporalNetworks","title":"TemporalNetworks.SpatTempMatching","text":"(:: SpatTempMatching)(::SpatTempMatching)(mlgraph :: MultilayerGraph{T}, norm :: Normalization, L_spat :: Matrix{Float64}, L_temp :: Matrix{Float64}) where T <: Multiplex\n\nComputes the diffusion constant a by matching the second spatial eigenvalue (first non-trivial spatial eigenvalue) to the first temporal eigenvalue for the supra-Laplacian corresponding to a multiplex temporal network. The multiplex structure guarantees a unique smallest intersection point. See [FroylandKoltai2023] and [AtnipFroylandKoltai2024] for details on this heuristic. \n\n\n\n\n\n","category":"type"},{"location":"#TemporalNetworks.SpatTempMatchingNonMultiplex","page":"TemporalNetworks","title":"TemporalNetworks.SpatTempMatchingNonMultiplex","text":"SpatTempMatchingNonMultiplex(evec_ind :: Int64, temp_ind :: Union{Int64, Nothng}, thresh :: Float64)\n\nCreates DiffusionEstimator instance to compute the diffusion constant a for the supra-Laplacian arising from a non-multiplex network by matching the eigenvalue corresponding to evec_ind to the corresponding multiplex temporal eigenvalue with index temp_ind. The parameter thresh is used by isspat(::MultilayerGraph{NonMultiplex}, ...) to check for spatial-like behaviour.\n\n(:: SpatTempMatchingNonMultiplex)(::SpatTempMatching)(mlgraph :: MultilayerGraph{T}, norm :: Normalization, L_spat :: Matrix{Float64}, L_temp :: Matrix{Float64}) where T <: NonMultiplex\n\nComputes the diffusion constant a using the matching heuristic for non-multiplex networks. See [FroylandKaliaKoltai2024] for details on this heuristic.\n\n\n\n\n\n","category":"type"},{"location":"#TemporalNetworks.bipartite_minimum_weight_edge_cover-Tuple{Any}","page":"TemporalNetworks","title":"TemporalNetworks.bipartite_minimum_weight_edge_cover","text":"Computes the minimum-weight edge cover of a bipartite graph G.  We assume there are m nodes on the left and n nodes on the right, where m≥n. The (i,j)th entry of the m×n matrix G contains the weight of the edge joining node i and node j. A binary m×n array is output, encoding the cover.\n\n\n\n\n\n","category":"method"},{"location":"#TemporalNetworks.bisection-Tuple{Any, Real, Real}","page":"TemporalNetworks","title":"TemporalNetworks.bisection","text":"bisection(f, a, b; fa = f(a), fb = f(b), ftol, wtol)\n\nBisection algorithm for finding the root f(x)  0 within the initial bracket [a,b].\n\nReturns a named tuple (x = x, fx = f(x), isroot = ::Bool, iter = ::Int, ismaxiter = ::Bool).\n\nTerminates when either\n\nabs(f(x)) < ftol (isroot = true),\nthe width of the bracket is ≤wtol (isroot = true, to account for discontinuities),\nmaxiter number of iterations is reached. (isroot = false, maxiter = true).\n\nwhich are tested for in the above order. Therefore, care should be taken not to make wtol too large.\n\n\n\n\n\n","category":"method"}]
}
