using Test
using CoinfectionSimulator
using Random
using DataFrames
using Distributions

# Set random seed for reproducible tests
Random.seed!(123)

@testset "CoinfectionSimulator.jl Tests" begin
    
    @testset "Basic Package Loading" begin
        @test isdefined(CoinfectionSimulator, :coinfection_simulator)
        @test isdefined(CoinfectionSimulator, :virtual_ecologist_sample)
        @test isdefined(CoinfectionSimulator, :prep_interaction_matrix)
    end
    
    include("test_simulator.jl")
    include("test_sampling.jl")
    include("test_data_prep.jl")
    include("test_utils.jl")
end
