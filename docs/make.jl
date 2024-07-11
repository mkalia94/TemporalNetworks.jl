import Pkg
Pkg.add("Documenter")
using Documenter
#push!(LOAD_PATH,"../")
push!(LOAD_PATH,"../src/")
using TemporalNetworks
makedocs(sitename="TemporalNetworks")
deploydocs(;repo = "github.com/mkalia94/TemporalNetworks.jl.git",target = "build",branch = "gh-pages")
