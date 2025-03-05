module TemporalNetworks

using LinearAlgebra, Plots, GraphRecipes, Base.Iterators, Statistics, Random, DataFrames, JuMP, HiGHS, PythonCall, CSV, StatsBase, SimpleWeightedGraphs, Parameters, ProgressMeter, CondaPkg, LaTeXStrings, StatsPlots, KrylovKit, SparseArrays 

include("struct/struct.jl")
include("struct/multilayer.jl")
include("struct/partition.jl")
include("struct/seba_partition.jl")
include("struct/block.jl")
include("utils/utils.jl")
include("utils/idl.jl")
include("utils/SEBA.jl")
include("utils/nonmultiplex.jl")
include("utils/plot.jl")
include("utils/compute_a.jl")
include("utils/bisection.jl")
include("gen/blocks.jl")
include("gen/blocks_nonmultiplex.jl")
include("utils/leiden/jc.jl")
include("utils/leiden/leiden.jl")
include("utils/leiden/BMEWCP.jl")

# Types
export Normalization, DegreeNormalization, IdentityNormalization
export TemporalConnectivity, Multiplex, NonMultiplex, NonMultiplexCompressed, NonMultiplexDense, DiffusionEstimator, SpatTempMatching, SpatTempMatchingNonMultiplex, RayleighBalancing
export MultilayerGraph

# Partitions
export SpectralPartition, lift
export SEBAPartition

# Block Graphs
export IntraBlockConnectivity, Block, BlockGraph, ScalingConnectivity, BlockGraphNonMultiplex 

# Utilities
export isspat, SEBA, cut, find_active, embed_compressed_evecs, plot, bisection, DynamicLaplacian, InflatedDynamicLaplacian

# Leiden
export leiden_slice, leiden_full

end # module TemporalNetworks
