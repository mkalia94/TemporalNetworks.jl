using Pkg
Pkg.activate("../.")
using TemporalNetworks
using Documenter, Example
makedocs(sitename="TemporalNetworks")
deploydocs(;repo = "github.com/mkalia94/TemporalNetworks.jl", devbranch="main")
