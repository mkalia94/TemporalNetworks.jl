using LinearAlgebra

function SEBA(V, Rinit=nothing)

# Inputs: 
# V is pxr matrix (r vectors of length p as columns)
# Rinit is an (optional) initial rotation matrix.

# Outputs:
# S is pxr matrix with columns approximately spanning the column space of V
# R is the optimal rotation that acts on V, which followed by thresholding, produces S

    maxiter = 5000   #maximum number of iterations allowed
    F = qr(V) # Enforce orthonormality
    V = Matrix(F.Q)
    p, r = size(V)
    μ = 0.99 / sqrt(p)

    S = zeros(size(V))
    # Perturb near-constant vectors
    for j = 1:r
        if maximum(V[:,j]) - minimum(V[:,j]) < 1e-14
            V[:,j] = V[:,j] .+ (rand(p, 1) .- 1 / 2) * 1e-12
        end
    end

    # Initialise rotation
    if Rinit ≡ nothing
        Rnew = Matrix(I, r, r)
    else
        # Ensure orthonormality of Rinit
        F = svd(Rinit)
        Rnew = F.U * F.Vt
    end

    R = zeros(r, r)
    iter = 0
    while norm(Rnew - R) > 1e-14 && iter < maxiter
        iter = iter + 1
        R = Rnew
        Z = V * R'
        # Threshold to solve sparse approximation problem
        for i = 1:r
            Si = sign.(Z[:,i]) .* max.(abs.(Z[:,i]) .- μ, zeros(p))
            S[:,i] = Si / norm(Si)
        end
        # Polar decomposition to solve Procrustes problem
        F = svd(S' * V, full=false)
        Rnew = F.U * F.Vt
    end

    # Choose correct parity of vectors and scale so largest value is 1
    for i = 1:r
        S[:,i] = S[:,i] * sign(sum(S[:,i]))
        S[:,i] = S[:,i] / maximum(S[:,i])
    end

    # Sort so that most reliable vectors appear first
    ind = sortperm(vec(minimum(S, dims=1)), rev=true)
    S = S[:, ind]

    return S, R

end


function double_basis_SEBA(graph :: MultilayerGraph{A}, evecs_for_SEBA :: Vector{Int64}, vecs :: Matrix{Float64}) where A
    new_basis = []

    for i in 1:length(evecs_for_SEBA)
        evec_ = reshape(vecs[:,evecs_for_SEBA[i]],graph.N,graph.T)
        norms = sqrt.(sum(abs2.(evec_),dims=1))
        new_vec = zeros(graph.N,graph.T)
        for j in 1:graph.T
            new_vec[:,j] = norms[j] * ones(graph.N)
        end
        push!(new_basis,new_vec)
    end
    new_basis = reshape.(new_basis,graph.N*graph.T,1) 
    basis = [vecs[:,evecs_for_SEBA] hcat(new_basis...)]
end


function _cut(partition :: SpectralPartition{MultilayerGraph{M},N}, V) where {M <: TemporalConnectivity, N <: Normalization}
    l = partition.graph.N * partition.graph.T
    V_comp = [x for x in Array(1:l) if x ∉ V]
    edges = collect(product(V, V_comp))
    Wfull = partition.L_spat + partition.a^2 .* partition.L_temp |> adj
    M <: Multiplex && return V_comp, edges, Wfull
    M <: NonMultiplex && return V_comp, edges, lift_adjacency(partition.graph, Wfull)
end

function lift_adjacency(graph :: MultilayerGraph{M}, W:: Matrix{Float64}) where {M <: NonMultiplex}
    W_lifted = zeros(graph.N*graph.T,graph.N*graph.T)
    ind_active = find_active(graph)
    ind_active = vcat([ind_active[i] .+ (i-1)*graph.N for i in 1:graph.T]...)
    W_lifted[ind_active, ind_active] = W
    W_lifted
end

function cut(partition :: SpectralPartition{MultilayerGraph{T}, IdentityNormalization}, V :: Vector) where T
    V_comp, edges, Wfull = _cut(partition, V)
    vol = minimum([length(V), length(V_comp)])
    return sum([Wfull[x,y] for x ∈ V, y ∈ V_comp ])/vol
end

function cut(partition :: SpectralPartition{MultilayerGraph{T}, DegreeNormalization}, V :: Vector) where T
    V_comp, edges, Wfull = _cut(partition, V)
    edges_V = sum(hcat([Wfull[x,:] for x in V]...))
    edges_V_comp = sum(hcat([Wfull[x,:] for x in V_comp]...))
    vol = minimum([edges_V, edges_V_comp])
    return sum([Wfull[x,y] for x in V, y in V_comp])/vol
end

function cuts_SEBA(partition :: SpectralPartition, vecs :: Matrix{Float64})
    
    cuts = zeros(size(vecs)[2])
    for i in 1:size(vecs)[2]
        V = findall(x->abs(x)>0.001,vecs[:,i])
        cuts[i] = cut(partition, V)
    end
    
    cuts
end

