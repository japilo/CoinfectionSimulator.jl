@testset "Virtual Ecologist Sampling Tests" begin
    
    @testset "Input validation tests" begin
        # Test empty virtual population
        @test_throws AssertionError virtual_ecologist_sample(
            virtual_population=Vector{Vector{Matrix{Bool}}}(),
            proportion_sampled=0.5,
            false_positive_rate=0.1,
            false_negative_rate=0.1
        )
        
        # Test invalid proportion_sampled
        valid_pop = [[rand(Bool, 10, 4) for _ in 1:5] for _ in 1:3]
        @test_throws AssertionError virtual_ecologist_sample(
            virtual_population=valid_pop,
            proportion_sampled=-0.1,
            false_positive_rate=0.1,
            false_negative_rate=0.1
        )
        
        @test_throws AssertionError virtual_ecologist_sample(
            virtual_population=valid_pop,
            proportion_sampled=1.1,
            false_positive_rate=0.1,
            false_negative_rate=0.1
        )
        
        # Test invalid false_positive_rate
        @test_throws AssertionError virtual_ecologist_sample(
            virtual_population=valid_pop,
            proportion_sampled=0.5,
            false_positive_rate=-0.1,
            false_negative_rate=0.1
        )
        
        # Test invalid false_negative_rate
        @test_throws AssertionError virtual_ecologist_sample(
            virtual_population=valid_pop,
            proportion_sampled=0.5,
            false_positive_rate=0.1,
            false_negative_rate=1.1
        )
    end
    
    @testset "Basic functionality tests" begin
        Random.seed!(42)
        
        # Create test data
        n_strains = 5
        n_individuals = 10
        n_timesteps = 3
        
        virtual_pop = Vector{Vector{Matrix{Bool}}}()
        for t in 1:n_timesteps
            timestep_data = Vector{Matrix{Bool}}()
            for i in 1:n_individuals
                individual_data = rand(Bool, n_strains, 4)
                push!(timestep_data, individual_data)
            end
            push!(virtual_pop, timestep_data)
        end
        
        # Test basic execution
        result = virtual_ecologist_sample(
            virtual_population=virtual_pop,
            proportion_sampled=0.5,
            false_positive_rate=0.0,
            false_negative_rate=0.0
        )
        
        @test size(result) == (n_timesteps, n_strains)
        @test eltype(result) == Bool
    end
    
    @testset "Edge case tests" begin
        Random.seed!(123)
        
        # Test with proportion_sampled = 0 (should still sample at least 1)
        virtual_pop = Vector{Vector{Matrix{Bool}}}([[rand(Bool, 5, 4) for _ in 1:10] for _ in 1:2])
        result = virtual_ecologist_sample(
            virtual_population=virtual_pop,
            proportion_sampled=0.0,
            false_positive_rate=0.0,
            false_negative_rate=0.0
        )
        @test size(result) == (2, 5)
        
        # Test with proportion_sampled = 1
        result = virtual_ecologist_sample(
            virtual_population=virtual_pop,
            proportion_sampled=1.0,
            false_positive_rate=0.0,
            false_negative_rate=0.0
        )
        @test size(result) == (2, 5)
        
        # Test with false_positive_rate = 1 and no true positives
        all_false_pop = Vector{Vector{Matrix{Bool}}}([[falses(3, 4) for _ in 1:5] for _ in 1:2])
        result = virtual_ecologist_sample(
            virtual_population=all_false_pop,
            proportion_sampled=1.0,
            false_positive_rate=1.0,
            false_negative_rate=0.0
        )
        @test all(result .== true)
        
        # Test with false_negative_rate = 1 and all true positives
        all_true_pop = [[trues(3, 4) for _ in 1:5] for _ in 1:2]
        # Set infected column (column 3) to true
        for timestep in all_true_pop
            for individual in timestep
                individual[:, 3] .= true
            end
        end
        result = virtual_ecologist_sample(
            virtual_population=all_true_pop,
            proportion_sampled=1.0,
            false_positive_rate=0.0,
            false_negative_rate=1.0
        )
        @test all(result .== false)
    end
    
    @testset "Output consistency tests" begin
        Random.seed!(456)
        
        # Test that output dimensions are consistent
        virtual_pop = [[rand(Bool, 8, 4) for _ in 1:12] for _ in 1:4]
        
        result1 = virtual_ecologist_sample(
            virtual_population=virtual_pop,
            proportion_sampled=0.3,
            false_positive_rate=0.05,
            false_negative_rate=0.1
        )
        
        result2 = virtual_ecologist_sample(
            virtual_population=virtual_pop,
            proportion_sampled=0.7,
            false_positive_rate=0.2,
            false_negative_rate=0.05
        )
        
        @test size(result1) == size(result2) == (4, 8)
        @test eltype(result1) == eltype(result2) == Bool
    end
end
