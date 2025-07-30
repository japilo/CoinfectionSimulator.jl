using Test
using CoinfectionSimulator
using Random
using DataFrames
using Distributions

# Set random seed for reproducible tests
Random.seed!(123)

@testset "CoinfectionSimulator.jl" begin

    @testset "Basic Package Loading" begin
        # Check core types are defined
        @test isdefined(CoinfectionSimulator, :DiseaseModel)
        @test isdefined(CoinfectionSimulator, :SIModel)
        @test isdefined(CoinfectionSimulator, :SIRModel)
        @test isdefined(CoinfectionSimulator, :SEIRModel)
        @test isdefined(CoinfectionSimulator, :SEIRSModel)
        @test isdefined(CoinfectionSimulator, :Individual)
        @test isdefined(CoinfectionSimulator, :Population)
        @test isdefined(CoinfectionSimulator, :SimulationParameters)
        @test isdefined(CoinfectionSimulator, :SamplingParameters)

        # Check core functions are defined
        @test isdefined(CoinfectionSimulator, :simulate)
        @test isdefined(CoinfectionSimulator, :sample)
        @test isdefined(CoinfectionSimulator, :create_interaction_matrix)

        # Check backward compatibility exports
        @test isdefined(CoinfectionSimulator, :coinfection_simulator)
        @test isdefined(CoinfectionSimulator, :virtual_ecologist_sample)
        @test isdefined(CoinfectionSimulator, :prep_interaction_matrix)
    end

    # Include test modules
    include("test_simulator.jl")
    include("test_sampling.jl")
    include("test_data_prep.jl")
    include("test_utils.jl")
end
