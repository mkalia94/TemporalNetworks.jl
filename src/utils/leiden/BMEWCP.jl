
"""Computes the minimum-weight edge cover of a bipartite graph `G`. 
We assume there are `m` nodes on the left and `n` nodes on the right, where `m≥n`.
The `(i,j)th` entry of the `m×n` matrix `G` contains the weight of the edge joining node `i` and node `j`.
A binary `m×n` array is output, encoding the cover."""
function bipartite_minimum_weight_edge_cover(G)

    #sample G below
    #G = [2 0; 2 0; 1 1]  Manus' example 12 Jan 2024
    #G = [2.5 2.5 2.5; 0 0 1; 1 1 1; 2 2 0]  Degenerate weights

    #ensure more rows than cols
    m, n = size(G)
    if m < n
        m, n = n, m
        G = G'
    end

    edgecover = Model(HiGHS.Optimizer)
    set_silent(edgecover)

    #define binary variable array
    @variable(edgecover, x[1:m, 1:n], Bin)

    # One cluster from the large group must be assigned to one cluster in the small group
    @constraint(edgecover, [i = 1:m], sum(x[i, :]) == 1)
    # Each cluster in the small group is assigned to at least one cluster in the large group
    @constraint(edgecover, [j = 1:n], sum(x[:, j]) ≥ 1)

    # Minimise total weight of the cover
    @objective(edgecover, Min, sum(G .* x))
    optimize!(edgecover)
    println("Objective value is ", objective_value(edgecover))

    return Bool.(value.(x))

end

