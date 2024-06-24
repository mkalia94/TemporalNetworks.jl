
function eigvecs_largemat(mat :: Matrix{Float64}, n :: Int64)
    dd = 2*maximum(mat)
    D = dd.*one(ones(size(mat)...))
    e_vals, e_vecs, info = eigsolve(D .- mat, rand(size(D)[2]), n, :LM)
    sum(abs2, info.normres) < 1e-10 ? (dd .- e_vals[1:n], hcat(e_vecs...)[:, 1:n]) : error("Krylov method didnt't converge.")
end


function (::SpatTempMatching)(mlgraph :: MultilayerGraph{T}, norm :: Normalization, L_spat :: Matrix{Float64}, L_temp :: Matrix{Float64}) where T <: Multiplex
    f = a -> (isspat(eigvecs(norm(L_spat + a^2 .* L_temp))[:,2], mlgraph.N, mlgraph.T, 0.01)) ? 1.0 : -1.0
    sol = bisection(f, 0.1, 50.0)
    !sol.ismaxiter ? sol.x : error("Max. iterations reached")
end


function (spnm::SpatTempMatchingNonMultiplex)(mlgraph :: MultilayerGraph{T}, norm :: Normalization, L_spat :: Matrix{Float64}, L_temp :: Matrix{Float64}) where T <: NonMultiplex
    @info "SpatTempMatching is not recommended for non-multiplex networks"
    f = a -> (abs(isspat(mlgraph, norm, L_spat, L_temp, a, spnm.evec_ind, spnm.temp_ind)) - spnm.thresh)
    sol = bisection(f, 0.1, 50.0)
    !sol.ismaxiter ? sol.x : error("Max. iterations reached")
end

#function (::SpatTempMatching)(mlgraph :: MultilayerGraph{T}, norm :: Normalization, L_spat :: Matrix{Float64}, L_temp :: Matrix{Float64}) where T <: NonMultiplex
#     error("SpatTempMatching doesn't work for non-multiplex graphs!")
#end

function _RayleighBalancing(mlgraph :: MultilayerGraph ,norm :: Normalization, L_spat :: Matrix{Float64}, L_temp :: Matrix{Float64}, a :: Float64)
    infL = L_spat + a^2 .* L_temp |> norm
    evecs = eigvecs(infL)
    rayleigh_spat = [sum((L_spat*evecs[:,i]) .* evecs[:,i]) for i in 1:10]
    rayleigh_temp  = [a^2 * sum((L_temp*evecs[:,i]) .* evecs[:,i]) for i in 1:10]
    rayleigh_spat .- rayleigh_temp
end


function (rb:: RayleighBalancing)(mlgraph :: MultilayerGraph ,norm :: Normalization, L_spat :: Matrix{Float64}, L_temp :: Matrix{Float64})
    f = a -> _RayleighBalancing(mlgraph, norm, L_spat, L_temp, a)[rb.ind]
    sol = bisection(f, 0.1, 50.0)
    if sol.isroot
        return sol.x
    else
        return 0.0
    end
end



