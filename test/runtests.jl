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
        @test isdefined(CoinfectionSimulator, :sample_populations)
        @test isdefined(CoinfectionSimulator, :create_interaction_matrix)



        # Check that internal helper functions are available
        @test isdefined(CoinfectionSimulator, :breeding!)
        @test isdefined(CoinfectionSimulator, :introduce_infections!)
        @test isdefined(CoinfectionSimulator, :add_individual!)
    end

    # Include test modules
    include("test_simulator.jl")
    include("test_sampling.jl")
    include("test_data_prep.jl")

end
