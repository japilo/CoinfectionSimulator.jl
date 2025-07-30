@testset "Simulator Tests" begin

    @testset "Types and Constructors" begin
        # Test disease model constructors
        @test SIModel(0.5, 0.1) isa SIModel
        @test SIRModel(0.5, 0.1, 0.2) isa SIRModel
        @test SEIRModel(0.5, 0.1, 0.2, 3) isa SEIRModel
        @test SEIRSModel(0.5, 0.1, 0.2, 3, 0.1) isa SEIRSModel

        # Test model validation
        @test_throws ArgumentError SIModel(-0.1, 0.1)  # Invalid transmission
        @test_throws ArgumentError SIModel(0.5, 1.1)   # Invalid mortality
        @test_throws ArgumentError SEIRModel(0.5, 0.1, 0.2, 0)  # Invalid latency

        # Test Individual constructor
        ind = Individual(BitMatrix([true false false false; true false false false]), 10)
        @test ind isa Individual
        @test ind.age == 10
        @test size(ind.state) == (2, 4)

        # Test Individual validation
        @test_throws ArgumentError Individual(BitMatrix([true false false]), 10)  # Wrong columns
        @test_throws ArgumentError Individual(BitMatrix([true true false false]), 10)  # Invalid state
        @test_throws ArgumentError Individual(BitMatrix([true false false false]), -1)  # Negative age

        # Test shorthand constructor
        ind2 = Individual(3, 5)
        @test ind2 isa Individual
        @test ind2.age == 5
        @test size(ind2.state) == (3, 4)
        @test all(ind2.state[:, 1])  # All susceptible
        @test !any(ind2.state[:, 2:4])  # Not in other states

        # Test Population constructor with different strain counts (should error)
        @test_throws ArgumentError Population([ind, ind2])

        pop = Population([ind, Individual(BitMatrix([true false false false; true false false false]), 20)])
        @test pop isa Population
        @test length(pop) == 2
    end

    @testset "Basic Simulation" begin
        # Create a small population for testing
        Random.seed!(123)
        ind1 = Individual(BitMatrix([false false true false; true false false false]), 10)
        ind2 = Individual(BitMatrix([true false false false; true false false false]), 10)
        pop = Population([ind1, ind2])

        # Simple SI model
        si_model = SIModel(0.8, 0.0)

        # Simulation parameters
        params = SimulationParameters(
            [si_model, si_model],   # Two strains with SI model
            [1.0 0.8; 0.8 1.0],     # Interaction matrix
            0.0,                    # No background mortality
            0.0,                    # No reproduction
            18,                     # Age of maturity
            :none,                  # No new introductions
            5                       # 5 time steps
        )

        # Run simulation
        results = simulate(pop, params)

        # Check results structure
        @test length(results) == 2  # Returns tuple of (populations, ages)
        @test length(results[1]) == 5  # 5 time steps
        @test length(results[2]) == 5  # 5 time steps

        # Check that infection spreads
        @test length(results[1][1]) == 2  # Initial population size
        @test count(ind -> ind[1, 3], results[1][1]) == 1  # One initially infected
        @test count(ind -> ind[1, 3], results[1][end]) >= 1  # At least one infected at end
    end

    @testset "Disease Model Dynamics" begin
        # Test SI model dynamics
        Random.seed!(456)
        pop_si = Population([Individual(BitMatrix([false false true false]), 20) for _ in 1:10])
        for _ in 1:90
            push!(pop_si.individuals, Individual(BitMatrix([true false false false]), 20))
        end

        params_si = SimulationParameters(
            [SIModel(0.3, 0.0)],  # High transmission, no mortality
            [1.0;;],              # No interaction (single strain)
            0.0,                  # No background mortality
            0.0,                  # No reproduction
            18,                   # Age of maturity
            :none,                # No new introductions
            10                    # 10 time steps
        )

        results_si = simulate(pop_si, params_si)

        # In SI model, infected stays infected
        final_infected = count(ind -> ind[1, 3], results_si[1][end])
        @test final_infected > 10  # Infection should spread

        # Test SIR model dynamics
        Random.seed!(457)
        pop_sir = Population([Individual(BitMatrix([false false true false]), 20) for _ in 1:10])
        for _ in 1:90
            push!(pop_sir.individuals, Individual(BitMatrix([true false false false]), 20))
        end

        params_sir = SimulationParameters(
            [SIRModel(0.3, 0.0, 0.2)],  # With recovery
            [1.0;;],
            0.0,
            0.0,
            18,
            :none,
            15
        )

        results_sir = simulate(pop_sir, params_sir)

        # In SIR model, we should see some recoveries
        final_recovered = count(ind -> ind[1, 4], results_sir[1][end])
        @test final_recovered > 0  # Some individuals should recover

        # Test SEIR model dynamics
        Random.seed!(458)
        pop_seir = Population([Individual(BitMatrix([false false true false]), 20) for _ in 1:10])
        for _ in 1:90
            push!(pop_seir.individuals, Individual(BitMatrix([true false false false]), 20))
        end

        params_seir = SimulationParameters(
            [SEIRModel(0.3, 0.0, 0.2, 3)],  # With latency
            [1.0;;],
            0.0,
            0.0,
            18,
            :none,
            20
        )

        results_seir = simulate(pop_seir, params_seir)

        # In SEIR model, we should see exposed individuals
        # Find maximum number of exposed at any time step
        max_exposed = maximum([count(ind -> ind[1, 2], pop) for pop in results_seir[1]])
        @test max_exposed > 0  # Should have some exposed individuals

        # Test SEIRS model dynamics
        Random.seed!(459)
        pop_seirs = Population([Individual(BitMatrix([false false true false]), 20) for _ in 1:10])
        for _ in 1:90
            push!(pop_seirs.individuals, Individual(BitMatrix([true false false false]), 20))
        end

        params_seirs = SimulationParameters(
            [SEIRSModel(0.3, 0.0, 0.2, 3, 0.1)],  # With immunity loss
            [1.0;;],
            0.0,
            0.0,
            18,
            :none,
            25
        )

        results_seirs = simulate(pop_seirs, params_seirs)

        # SEIRS model should have both recovered and newly susceptible individuals
        final_recovered = count(ind -> ind[1, 4], results_seirs[1][end])
        @test final_recovered > 0  # Some should be recovered

        # Note: Testing for immunity loss is challenging with random processes
        # We'd need a more controlled test to verify that specific behavior
    end

    @testset "Demographic Processes" begin
        # Test mortality
        Random.seed!(567)
        pop = Population([Individual(BitMatrix([false false true false]), 20) for _ in 1:100])

        params = SimulationParameters(
            [SIModel(0.0, 0.2)],  # No transmission, high disease mortality
            [1.0;;],
            0.0,                  # No background mortality
            0.0,                  # No reproduction
            18,
            :none,
            10
        )

        results = simulate(pop, params)

        # Population should decrease due to mortality
        @test length(results[1][end]) < 100

        # Test reproduction
        Random.seed!(568)
        pop = Population([Individual(BitMatrix([true false false false]), 25) for _ in 1:50])

        params = SimulationParameters(
            [SIModel(0.0, 0.0)],  # No disease dynamics
            [1.0;;],
            0.0,                  # No mortality
            0.5,                  # High fecundity
            18,                   # Mature at age 18
            :none,
            10
        )

        results = simulate(pop, params)

        # Population should increase due to reproduction
        @test length(results[1][end]) > 50

        # Newborns should be age 0 or small age in the final population
        @test any(ind -> ind.age < 10, results[1][end])
    end

    @testset "Multi-strain Interactions" begin
        # Test that interactions affect transmission when both strains are present
        Random.seed!(12345)

        # Create a population with some individuals infected with each strain
        pop = Population(Individual[])

        # Add individuals infected with strain 1 only
        for _ in 1:40
            push!(pop.individuals, Individual(BitMatrix([false false true false; true false false false]), 20))
        end

        # Add a few individuals infected with strain 2 only (to seed strain 2 transmission)
        for _ in 1:5
            push!(pop.individuals, Individual(BitMatrix([true false false false; false false true false]), 20))
        end

        # Add susceptible individuals
        for _ in 1:55
            push!(pop.individuals, Individual(BitMatrix([true false false false; true false false false]), 20))
        end

        # Test with no facilitation first
        interactions_baseline = [1.0 1.0; 1.0 1.0]
        params_baseline = SimulationParameters(
            [SIModel(0.1, 0.0), SIModel(0.1, 0.0)],
            interactions_baseline,
            0.0, 0.0, 18, :none, 30
        )

        results_baseline = simulate(pop, params_baseline)
        strain2_baseline = count(ind -> ind[2, 3], results_baseline[1][end])

        # Now test with facilitation (strain 1 facilitates strain 2 transmission)
        Random.seed!(12345)  # Reset seed for comparison

        # Create identical population
        pop_facilitated = Population(Individual[])
        for _ in 1:40
            push!(pop_facilitated.individuals, Individual(BitMatrix([false false true false; true false false false]), 20))
        end
        for _ in 1:5
            push!(pop_facilitated.individuals, Individual(BitMatrix([true false false false; false false true false]), 20))
        end
        for _ in 1:55
            push!(pop_facilitated.individuals, Individual(BitMatrix([true false false false; true false false false]), 20))
        end

        interactions_facilitated = [1.0 1.0; 2.0 1.0]  # Strain 1 facilitates strain 2
        params_facilitated = SimulationParameters(
            [SIModel(0.1, 0.0), SIModel(0.1, 0.0)],
            interactions_facilitated,
            0.0, 0.0, 18, :none, 30
        )

        results_facilitated = simulate(pop_facilitated, params_facilitated)
        strain2_facilitated = count(ind -> ind[2, 3], results_facilitated[1][end])

        # With facilitation, strain 2 should spread more (or at least as much)
        @test strain2_facilitated >= strain2_baseline

        # Simple test: verify that interactions don't break the basic simulation
        @test length(results_facilitated[1][end]) > 0  # Population should still exist
        @test sum(count(ind -> ind[i, j], results_facilitated[1][end]) for i in 1:2, j in 1:4) > 0  # Some disease states should exist
    end

    @testset "Backward Compatibility" begin
        # Test that the old API still works
        # Create an identical setup for both APIs with a fixed seed

        # Create identical initial conditions
        init_pop_legacy = [Matrix{Bool}([true false false false; true false false false]) for _ in 1:20]
        init_pop_legacy[1][1, 1] = false  # First individual, strain 1 infected
        init_pop_legacy[1][1, 3] = true

        # Create fixed ages and parameters
        init_ages = ones(Int, 20)

        # Use a single infected individual to start, no randomness in the model
        Random.seed!(12345)
        results_legacy = coinfection_simulator(
            initial_pop=init_pop_legacy,
            ages=init_ages,
            interactions=[1.0 0.8; 0.8 1.0],
            disease_type=["si", "si"],
            base_mortality=0.0,     # No mortality
            disease_mortality=[0.0, 0.0],
            fecundity=0.0,          # No births
            transmission=[0.0, 0.0], # No transmission - should preserve initial infected count
            time_steps=5,
            age_maturity=1,
            introduction="none"     # Explicit setting
        )

        # Check that the results have the expected structure
        @test length(results_legacy) == 2
        @test length(results_legacy[1]) == 5
        @test length(results_legacy[2]) == 5

        # Create equivalent setup for new API
        Random.seed!(12345)

        # Build identical population
        ind1 = Individual(BitMatrix([false false true false; true false false false]), 1)
        pop_new = Population([ind1])
        for _ in 1:19
            push!(pop_new.individuals, Individual(BitMatrix([true false false false; true false false false]), 1))
        end

        params_new = SimulationParameters(
            [SIModel(0.0, 0.0), SIModel(0.0, 0.0)],  # No transmission
            [1.0 0.8; 0.8 1.0],
            0.0,                    # No mortality
            0.0,                    # No births
            1,
            :none,                  # No introductions
            5
        )

        results_new = simulate(pop_new, params_new)

        # Compare results
        @test length(results_new[1]) == length(results_legacy[1])
        @test all(results_new[2][i] == results_legacy[2][i] for i in 1:5)

        # Check that the final infection counts match - should remain at 1 since no transmission
        infected_legacy = sum(sum(m[1, 3] for m in results_legacy[1][end]))
        infected_new = count(ind -> ind[1, 3], results_new[1][end])
        @test infected_legacy == infected_new
        @test infected_legacy == 1  # Only the initial infection
    end
end
