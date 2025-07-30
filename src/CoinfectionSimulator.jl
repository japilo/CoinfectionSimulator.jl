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
include("simulator.jl")
include("sampling.jl")
include("data_prep.jl")

# Export main types
export DiseaseModel, SIModel, SIRModel, SEIRModel, SEIRSModel
export Population, Individual, SimulationParameters, SamplingParameters

# Export main functions
export simulate, sample_populations, create_interaction_matrix

# Convenience re-exports of types and functions with legacy names for backward compatibility
export coinfection_simulator, virtual_ecologist_sample, prep_interaction_matrix

end # module CoinfectionSimulator
