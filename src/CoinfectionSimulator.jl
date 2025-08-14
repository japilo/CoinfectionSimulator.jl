module CoinfectionSimulator

# Standard library imports
using Random
using LinearAlgebra

# External dependencies
using Distributions
using StatsBase
using DataFrames

# Type definitions
include("types.jl")

# Core functionality
include("simulation_helpers.jl")
include("simulator.jl")
include("sampling.jl")
include("data_prep.jl")

# Export main types
export DiseaseModel, SIModel, SIRModel, SEIRModel, SEIRSModel
export Population, Individual, SimulationParameters, SamplingParameters

# Export main functions
export simulate, sample_populations, create_interaction_matrix



end # module CoinfectionSimulator
