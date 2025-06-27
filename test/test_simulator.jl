@testset "Coinfection Simulator Tests" begin

    @testset "Input validation tests" begin
        # Test empty population
        @test_throws AssertionError coinfection_simulator(
            initial_pop=Vector{Matrix{Bool}}(),
            ages=Vector{Int}(),
            interactions=Matrix{Float64}(undef, 0, 0),
            disease_type=String[],
            base_mortality=0.0,
            disease_mortality=Float64[],
            fecundity=0.0,
            transmission=Float64[],
            time_steps=1,
            age_maturity=1
        )
        
        # Test Matrix{Bool} input
        pop = [Matrix{Bool}([true false false false; true false false false]) for _ in 1:10]
        pop[1][1, 1] = false
        pop[1][1, 3] = true
        
        result = coinfection_simulator(
            initial_pop=pop,
            ages=ones(Int, 10),
            interactions=[1.0 0.9; 0.9 1.0],
            disease_type=["si", "si"],
            base_mortality=0.0,
            disease_mortality=[0.0, 0.0],
            fecundity=0.0,
            transmission=[0.5, 0.5],
            time_steps=5,
            age_maturity=1
        )
        
        @test length(result) == 2  # Returns tuple of (populations, ages)
        @test length(result[1]) == 5  # 5 time steps
        @test length(result[2]) == 5  # 5 time steps

        # Test basic valid input with BitMatrix
        pop_bit = [BitMatrix([true false false false; true false false false]) for _ in 1:10]
        pop_bit[1][1, 1] = false
        pop_bit[1][1, 3] = true
        
        result_bit = coinfection_simulator(
            initial_pop=pop_bit,
            ages=ones(Int, 10),
            interactions=[1.0 0.9; 0.9 1.0],
            disease_type=["si", "si"],
            base_mortality=0.0,
            disease_mortality=[0.0, 0.0],
            fecundity=0.0,
            transmission=[0.5, 0.5],
            time_steps=5,
            age_maturity=1
        )
        
        @test length(result_bit) == 2  # Returns tuple of (populations, ages)
        @test length(result_bit[1]) == 5  # 5 time steps
        @test length(result_bit[2]) == 5  # 5 time steps
    end

    @testset "BitMatrix compatibility tests" begin
        Random.seed!(789)
        
        # Test that BitMatrix and Matrix{Bool} produce equivalent results
        pop_bool = [Matrix{Bool}([true false false false; true false false false]) for _ in 1:20]
        pop_bool[1][1, 1] = false
        pop_bool[1][1, 3] = true
        
        pop_bit = [BitMatrix([true false false false; true false false false]) for _ in 1:20]
        pop_bit[1][1, 1] = false
        pop_bit[1][1, 3] = true
        
        # Run simulation with identical parameters
        Random.seed!(123)
        result_bool = coinfection_simulator(
            initial_pop=pop_bool,
            ages=ones(Int, 20),
            interactions=[1.0 0.8; 0.8 1.0],
            disease_type=["si", "si"],
            base_mortality=0.0,
            disease_mortality=[0.0, 0.0],
            fecundity=0.0,
            transmission=[0.3, 0.3],
            time_steps=5,
            age_maturity=1
        )
        
        Random.seed!(123)
        result_bit = coinfection_simulator(
            initial_pop=pop_bit,
            ages=ones(Int, 20),
            interactions=[1.0 0.8; 0.8 1.0],
            disease_type=["si", "si"],
            base_mortality=0.0,
            disease_mortality=[0.0, 0.0],
            fecundity=0.0,
            transmission=[0.3, 0.3],
            time_steps=5,
            age_maturity=1
        )
        
        # Results should be equivalent
        @test length(result_bool[1]) == length(result_bit[1])
        @test result_bool[2] == result_bit[2]  # Ages should be identical
        
        # Check that population states are equivalent
        for t in 1:length(result_bool[1])
            @test length(result_bool[1][t]) == length(result_bit[1][t])
            for i in 1:length(result_bool[1][t])
                @test Matrix{Bool}(result_bool[1][t][i]) == Matrix{Bool}(result_bit[1][t][i])
            end
        end
    end

    @testset "Disease model validation" begin
        pop = [Matrix{Bool}([true false false false;]) for _ in 1:10]
        
        # Test SIR model requires recovery parameter
        @test_throws AssertionError coinfection_simulator(
            initial_pop=pop,
            ages=ones(Int, 10),
            interactions=[1.0;;],
            disease_type=["sir"],
            base_mortality=0.0,
            disease_mortality=[0.0],
            fecundity=0.0,
            transmission=[0.5],
            time_steps=2,
            age_maturity=1,
            # Missing recovery parameter
        )
        
        # Test SEIR model requires latency parameter
        @test_throws AssertionError coinfection_simulator(
            initial_pop=pop,
            ages=ones(Int, 10),
            interactions=[1.0;;],
            disease_type=["seir"],
            base_mortality=0.0,
            disease_mortality=[0.0],
            fecundity=0.0,
            transmission=[0.5],
            time_steps=2,
            age_maturity=1,
            recovery=[0.1]
            # Missing latency parameter
        )
        
        # Test SEIRS model requires immunity_loss parameter
        @test_throws AssertionError coinfection_simulator(
            initial_pop=pop,
            ages=ones(Int, 10),
            interactions=[1.0;;],
            disease_type=["seirs"],
            base_mortality=0.0,
            disease_mortality=[0.0],
            fecundity=0.0,
            transmission=[0.5],
            time_steps=2,
            age_maturity=1,
            recovery=[0.1],
            latency=[5]
            # Missing immunity_loss parameter
        )
    end

    @testset "Disease progression tests" begin
        Random.seed!(456)
        
        # Test SI model
        pop = [Matrix{Bool}([true false false false; true false false false]) for _ in 1:100]
        pop[1][1, 1] = false
        pop[1][1, 3] = true  # Start with one infected individual
        
        results = coinfection_simulator(
            initial_pop=pop,
            ages=ones(Int, 100),
            interactions=[1.0 0.5; 0.5 1.0],
            disease_type=["si", "si"],
            base_mortality=0.0,
            disease_mortality=[0.0, 0.0],
            fecundity=0.0,
            transmission=[0.8, 0.8],
            time_steps=10,
            age_maturity=1
        )
        
        # Check that infection can spread
        initial_infected = sum(m[1,3] for m in results[1][1])
        final_infected = sum(m[1,3] for m in results[1][end])
        @test final_infected >= initial_infected
    end

    @testset "Edge cases" begin
        # Test with mortality = 0 and fecundity = 0 (stable population)
        pop = [Matrix{Bool}([true false false false;]) for _ in 1:5]
        
        results = coinfection_simulator(
            initial_pop=pop,
            ages=ones(Int, 5),
            interactions=[1.0;;],
            disease_type=["si"],
            base_mortality=0.0,
            disease_mortality=[0.0],
            fecundity=0.0,
            transmission=[0.0],  # No transmission
            time_steps=3,
            age_maturity=1,
            introduction="none"
        )
        
        # Population should remain stable
        @test length(results[1][1]) == length(results[1][end])
        @test all(results[2][end] .== results[2][1] .+ 2)  # Ages should increase by 2
    end
end
