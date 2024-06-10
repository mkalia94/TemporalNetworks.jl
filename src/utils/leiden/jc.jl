"""
This program is to discover and plot the communities of a network by leiden algorithm
NOTE: the leiden algorithm is implemented by the python package leidenalg, so
        before doing community discovery, Conda, PyCall have to be installed as
        follows:

            import Pkg
            Pkg.add("Conda")
            Pkg.add("PyCall")
            Pkg.build("PyCall")
            
            using Conda

            Conda.pip_interop(true)
            #Conda.pip("install", "scipy")
            #Conda.pip("install", "numpy")
            Conda.pip("install", "leidenalg")

Contributors: Xiaoshan Nian
Date: August, 2020
Email: cen@njust.edu.cn
Github: https://github.com/yottoo/JuliaCommunity
"""

module JuliaCommunity

using PyCall
using DataFrames, CSV
using Statistics, StatsBase, Random
using Parameters, ProgressMeter, SimpleWeightedGraphs

export JuliaCommunityInstance, 
        discover_communities, 
        optimise_resolution, 
        plot_community

const leiden = pyimport("leidenalg")
const ig = pyimport("igraph")

@with_kw mutable struct PartitionMethod
    louvain::String = "louvain"
    CPM::String = "CPM"
    modularity::String = "modularity"
end


@with_kw mutable struct JuliaCommunityInstance
    task_series::String = ""

    methods::PartitionMethod = PartitionMethod()
    method::String = "CPM"    # CPM, Modularity (for leiden algorithm) and Louvain
    """
    NETWORK: the network data.
        from: id of the node from
        to: id of the node to
        weight: edge weight (if edge_weighted is false, this column could be ignored)
    NODES: the data for nodes, with ID required
        id: node id
    NODE_IMPORTANCES: the importances of vertices
        id: node id
        importance
    """    
    network::DataFrame = DataFrame()
    graph = nothing
    igraph = nothing
    nodes::DataFrame = DataFrame()
    max_node_id::Int = 1
    node_importance_field::String = "importance"
    node_label_field::String = "id"
    
    edge_weighted::Bool = true
    node_weighted::Bool = false
    is_directed::Bool = true
    γ::Float16 = 0.001    

    """
    COMMUNITIES: communities discovered.
        c: community id (start from 1)
        size: community size
        cluster_cof: clustering coefficent of the community
    MEMBERSHIPS: memberships indicating the nodes belonging to communities
        id: node id
        c: community id   
    """
    communities::DataFrame = DataFrame()
    membership_vector::Array{Int} = []
    memberships::DataFrame = DataFrame()
    n_community::Int = 0

    modularity::Float64 = 0
    quality::Float64 = 0
end


function JuliaCommunityInstance(network::DataFrame; 
        nodes::DataFrame=DataFrame(), 
        node_label_field::String="id",
        node_importance_field="importance", 
        edge_weighted::Bool=true, 
        node_weighted::Bool=false, 
        is_directed::Bool=true, 
        method::String="CPM", 
        task_series::String="",)
    filter!(row -> row.from > 0 && row.to > 0, network)    
    jc = JuliaCommunityInstance()
    jc.task_series = replace("_" * replace(task_series, "-" => "_"), "__" => "_")
    set_method(jc, method=method)
    #print("\nLoading the network and nodes data...")
    jc.network = network
    check_network(jc)
    jc.nodes = nodes
    check_nodes(jc)
    jc.max_node_id = maximum(jc.nodes.id)
    jc.node_importance_field = node_importance_field
    jc.node_label_field = node_label_field
    jc.edge_weighted = edge_weighted
    jc.node_weighted = node_weighted
    jc.is_directed = is_directed
    
    if edge_weighted
        filter!(row -> row.weight > 0, jc.network) 
    #= ============================================
    else if !is_directed
        jc.network = vcat(jc.network, jc.network)
        unique(jc.network)
    ============================================ =#
    end

    #print("\n\tBuilding the graph from the network...")
    if jc.is_directed && jc.edge_weighted
        jc.graph = SimpleWeightedDiGraph(jc.network.from, jc.network.to, jc.network.weight)
        jc.igraph = ig.Graph(zip(jc.network.from .- 1, jc.network.to .- 1), directed=true, edge_attrs=Dict("weight" => jc.network.weight))
    elseif jc.is_directed        
        jc.graph = SimpleDiGraph(jc.max_node_id)
        for i in 1:nrow(jc.network) add_edge!(jc.graph, jc.network.from[i], jc.network.to[i]) end
        jc.igraph = ig.Graph(zip(jc.network.from .- 1, jc.network.to .- 1), directed=true)
    elseif jc.edge_weighted
        jc.graph = SimpleWeightedGraph(jc.network.from, jc.network.to, jc.network.weight)
        jc.igraph = ig.Graph(zip(jc.network.from .- 1, jc.network.to .- 1), directed=false, edge_attrs=Dict("weight" => jc.network.weight))
    else
        jc.graph = SimpleGraph(jc.max_node_id)
        for i in 1:nrow(jc.network) add_edge!(jc.graph, jc.network.from[i], jc.network.to[i]) end
        jc.igraph = ig.Graph(zip(jc.network.from .- 1, jc.network.to .- 1), directed=false)
    end
      
    jc        
end


function check_network(jc::JuliaCommunityInstance)
    columns = names(jc.network)
    if isnothing(findfirst(name -> name == "from", columns)) throw(error("The network data must have a 'from' column.")) end
    if isnothing(findfirst(name -> name == "to", columns)) throw(error("The network data must have a 'to' column.")) end
    if jc.edge_weighted && isnothing(findfirst(name -> name == "weight", columns)) throw(error("The network data must have a 'weight' column.")) end    
end


function check_nodes(jc::JuliaCommunityInstance)
    columns = names(jc.nodes)
    if isnothing(findfirst(name -> name == "id", columns)) throw(error("The nodes data must have a 'id' column.")) end    
end

function set_method(jc::JuliaCommunityInstance; method::String="CPM")
    jc.method = method
    if (jc.method != jc.methods.louvain && jc.method != jc.methods.CPM && jc.method != jc.methods.modularity) 
        throw(error("The partition method has to be louvain, CMP or modularity.")) 
    end
end

"""
_louvain:    the louvain algorithm implemented by leiden package    
"""
function _louvain(g)
    optimiser = leiden.Optimiser()
    partitions = leiden.ModularityVertexPartition(g)
    partitions_agg = partitions.aggregate_partition()
    while optimiser.move_nodes(partitions) > 0
        partitions.from_coarse_partition(partitions_agg)
        partitions_agg = partitions_agg.aggregate_partition()
    end
    partitions
end


"""
discover_communities:    discover the communities by leiden algorithm (both CMP and modularity method) 
    and louvain algorithm.
"""
function discover_communities(jc::JuliaCommunityInstance; mute::Bool=true)
    if !mute print("\nDiscovering the communities for the built network......") end
    
    partitions = nothing
    if jc.method == jc.methods.CPM
        partitions = leiden.find_partition(jc.igraph, leiden.CPMVertexPartition, resolution_parameter=jc.γ)
    elseif jc.method == jc.methods.modularity
        partitions = leiden.find_partition(jc.igraph, leiden.ModularityVertexPartition)
    elseif jc.method == jc.methods.louvain
        partitions = _louvain(jc.igraph)
    end
    
    if isnothing(partitions) return end

    jc.n_community = length(partitions)
    if !mute println("\t\t$(jc.n_community) communities have been discovered.") end
    
    # partitions = leiden.find_partition(g, leiden.ModularityVertexPartition, resolution_parameter=0.2)
    # The following code could not be done to Leiden community members 
    # for i in 1:length(partitions) partitions[i] .+= 1 end
    # partitions.membership .+= 1
    jc.membership_vector = partitions.membership .+ 1
    # println(partitions.membership)
    
    #jc.modularity = modularity(jc.graph, jc.membership_vector)
    jc.quality =  jc.method == jc.methods.CPM ? partitions.quality() : jc.modularity
    
    jc.memberships = DataFrame(id=1:maximum(vcat(jc.network.from, jc.network.to)), c=jc.membership_vector)
    # jc.memberships = DataFrame(id=1:nv(g), c=jc.membership_vector)
    jc.communities = DataFrame(c=1:jc.n_community, size=length.(partitions))
end

function optimise_resolution(jc::JuliaCommunityInstance; γ_from::Float64=0.0001, γ_end::Float64=0.01, γ_step::Float64=0.0001)
    println("\n")
    qualities = DataFrame(resolution=[], n_community=[], modularity=[], quality=[])
    jc_copy = deepcopy(jc)
    progress = Progress(length(γ_from:γ_step:γ_end), desc="Finding the best resolution γ for the Leiden-based community discovery algorithm: ")
    for γ  in γ_from:γ_step:γ_end
        jc_copy.γ = γ
        discover_communities(jc_copy, mute = true)
        push!(qualities, (γ, jc_copy.n_community, jc_copy.modularity, jc_copy.quality))
        next!(progress)
        #print("\n\t\tResolution: $γ: $(jc_copy.n_community) commxunities discovered; Modularity: $(jc_copy.modularity); CPM Quality: $(jc_copy.quality).") 
    end
    CSV.write("data/community_discover_optimisation-$(jc.method)$(jc.task_series).csv", qualities)

    p_modularity = plot(
        layer(qualities, x=:resolution, y=:modularity, Geom.line, Geom.point, Theme(default_color="blue")),        
        Guide.xticks(ticks=γ_from:γ_step * 2:γ_end),
        Guide.xlabel("resolution γ"),
        Guide.ylabel("modularity"),
        Theme(major_label_font_size=10pt),
        Scale.y_continuous(format=:plain)
    )
    fig_modularity = "fig/community_discover_optimisation-modularities-$(jc.method)$(jc.task_series).svg"
    draw(SVG(fig_modularity, 24cm, 16cm), p_modularity);
    open_file(fig_modularity)

    p_quality = plot(
        layer(qualities, x=:resolution, y=:quality, Geom.line, Geom.point, Theme(default_color="blue")),        
        Guide.xticks(ticks=γ_from:γ_step * 2:γ_end),
        Guide.xlabel("resolution γ"),
        Guide.ylabel("quality"),
        Theme(major_label_font_size=10pt),
        Scale.y_continuous(format=:plain)
    )
    fig_quality = "fig/community_discover_optimisation-qualities-$(jc.method)$(jc.task_series).svg"
    draw(SVG(fig_quality, 24cm, 16cm), p_quality);
    open_file(fig_quality)

    p_n_communities = plot(
        layer(qualities, x=:resolution, y=:n_community, Geom.line, Geom.point, Theme(default_color="blue")),        
        Guide.xticks(ticks=γ_from:γ_step * 2:γ_end),
        Guide.xlabel("resolution γ"),
        Guide.ylabel("number of communities"),
        Theme(major_label_font_size=10pt),
        Scale.y_continuous(format=:plain)
    )
    fig_n_communities = "fig/community_discover_optimisation-n-communities-$(jc.method)$(jc.task_series).svg"
    draw(SVG(fig_n_communities, 24cm, 16cm), p_n_communities);
    open_file(fig_n_communities)
end


end
