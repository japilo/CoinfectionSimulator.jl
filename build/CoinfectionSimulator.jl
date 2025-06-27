module CoinfectionSimulator

using Random
using Distributions
using StatsBase
using DataFrames

# Export main functions
export coinfection_simulator, virtual_ecologist_sample, prep_interaction_matrix

# Include submodules
include("simulator.jl")
include("sampling.jl")
include("data_prep.jl")
include("utils.jl")

end # module CoinfectionSimulator
