function directsum(vec :: Vector) # vector of matrices
    M = length(vec)
    N = [size(vec[i])[1] for i in 1:M]
    W = zeros(sum(N),sum(N))

    ctr = 0
    for i in 1:M
        W[ctr+1:ctr + N[i],ctr+1:ctr + N[i]] = vec[i]
        ctr = ctr + N[i]
    end
    W
end

function directsum(vec :: Vector{SparseMatrixCSC})
    M = length(vec)
    N = [vec[i].n for i in 1:M]
    W = spzeros(sum(N),sum(N))
    ctr = 0
    for i in 1:M
        W[ctr+1:ctr + N[i],ctr+1:ctr + N[i]] = vec[i]
        ctr = ctr + N[i]
    end
    W
end



function lap(W :: Matrix{Float64})
    return diagm(sum(abs.(W),dims=1)[1,:]) - W
end

function lap(W :: SparseMatrixCSC)
    return spdiagm(sum(abs.(W), dims=1)[1,:]) - W
end

function adj(L :: Matrix{Float64})
    diagm(diag(L)) .- L
end

function adj(L :: SparseMatrixCSC)
    spdiagm(diag(L)) .- L
end

function edgelist(W :: Union{Matrix, SparseMatrixCSC})
    edges = []
    for i in 1: size(W)[1], j in 1:size(W)[2]
        if W[i,j] > 0
            push!(edges,[i,j])
        end
    end
    edges
end

function isspat(evec :: Vector{Float64}, N :: Int64, T :: Int64, thresh :: Float64)
    vec = reshape(evec, N, T)
    vec = vec ./ sqrt.(sum(abs2, vec, dims=1))
    flag = var(sum(vec, dims=1))
    flag < thresh
end


function color_vector(vec)
    vec_colours = [colorant"rgba(0,114,178,1)",colorant"rgba(213,94,0,1)",colorant"rgba(231,159,0,1)",colorant"rgba(1,158,115,1)",colorant"rgba(159,95,159,1)",colorant"rgba(240,228,65,1)", colorant"rgba(87,180,233,1)",colorant"rgba(237,30,36,1)",colorant"rgba(3,255,0,1)",colorant"rgba(196,0,255,1)",colorant"rgba(125,125,125,1)"]
    if 0 ∉ vec
        return vec_colours[Int64.(vec)]
    else
        vec_temp = [colorant"rgba(0,114,178,1)" for i in 1:length(vec)] 
        ind = findall(x->x ≠ 0, vec)
        if !isempty(ind)
            vec_temp[ind] = vec_colours[Int64.(vec[ind])]
        end
        ind_others = [x for x in 1:length(vec) if x ∉ ind]
        vec_temp[ind_others] .= colorant"rgba(255,255,255,1)"
        return [vec_temp...]
    end
        
end

function get_ticks(M :: Int64)

    x1 = 0
    xn = M
    x2 = floor((xn-x1)/4) |> Int64
    x3 = 2*x2 |> Int64
    x4 = 3*x2 |> Int64
    x5 = 4*x2 |> Int64
    return [x2,x3,x4,x5]
end

