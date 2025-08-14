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
        @test length(results) == 5  # 5 time steps

        # Check that infection spreads
        @test length(results[1]) == 2  # Initial population size
        @test count(ind -> ind[1, 3], results[1]) == 1  # One initially infected
        @test count(ind -> ind[1, 3], results[end]) >= 1  # At least one infected at end
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
        final_infected = count(ind -> ind[1, 3], results_si[end])
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
        final_recovered = count(ind -> ind[1, 4], results_sir[end])
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
        max_exposed = maximum([count(ind -> ind[1, 2], pop) for pop in results_seir])
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
        final_recovered = count(ind -> ind[1, 4], results_seirs[end])
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
        @test length(results[end]) < 100

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
        @test length(results[end]) > 50

        # Newborns should be age 0 or small age in the final population
        @test any(ind -> ind.age < 10, results[end])
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
        strain2_baseline = count(ind -> ind[2, 3], results_baseline[end])

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
        strain2_facilitated = count(ind -> ind[2, 3], results_facilitated[end])

        # With facilitation, strain 2 should spread more (or at least as much)
        @test strain2_facilitated >= strain2_baseline

        # Simple test: verify that interactions don't break the basic simulation
        @test length(results_facilitated[end]) > 0  # Population should still exist
        @test sum(count(ind -> ind[i, j], results_facilitated[end]) for i in 1:2, j in 1:4) > 0  # Some disease states should exist
    end

    @testset "Base Mortality Applied Only Once" begin
        # Test that base mortality is applied only once per individual per time step,
        # regardless of the number of strains
        Random.seed!(12345)

        # Create a large population to reduce stochastic variation
        n_individuals = 1000
        n_strains = 3  # Multiple strains to test the fix
        base_mortality_rate = 0.05  # 5% base mortality per time step
        n_time_steps = 20

        # Create population - all individuals start susceptible
        pop = Population([Individual(n_strains, 20) for _ in 1:n_individuals])

        # Set up simulation with NO disease mortality, NO transmission, NO reproduction
        models = [SIModel(0.0, 0.0) for _ in 1:n_strains]  # No transmission, no disease mortality

        params = SimulationParameters(
            models,
            ones(Float64, n_strains, n_strains),  # No strain interactions
            base_mortality_rate,     # Only base mortality
            0.0,                     # No reproduction
            100,                     # High maturity age (no reproduction)
            :none,                   # No strain introductions
            n_time_steps
        )

        # Run simulation
        results = simulate(pop, params)
        populations = results

        # Calculate actual mortality rates
        initial_pop_size = length(populations[1])
        final_pop_size = length(populations[end])

        # Expected survival probability over n_time_steps: (1 - base_mortality)^n_time_steps
        expected_survival_prob = (1 - base_mortality_rate)^(n_time_steps - 1)
        expected_final_size = initial_pop_size * expected_survival_prob

        # Calculate tolerance (allow for 2 standard deviations of binomial distribution)
        # For binomial with n trials and probability p, std = sqrt(n * p * (1-p))
        std_dev = sqrt(initial_pop_size * (1 - expected_survival_prob) * expected_survival_prob)
        tolerance = 2 * std_dev

        # Test that final population size is within expected range
        @test abs(final_pop_size - expected_final_size) <= tolerance

        # Also test that we actually had some deaths (mortality > 0)
        @test final_pop_size < initial_pop_size

        # Test with a single strain to make sure we get the same result
        Random.seed!(12345)  # Same seed for comparison

        pop_single = Population([Individual(1, 20) for _ in 1:n_individuals])
        models_single = [SIModel(0.0, 0.0)]

        params_single = SimulationParameters(
            models_single,
            ones(Float64, 1, 1),
            base_mortality_rate,
            0.0,
            100,
            :none,
            n_time_steps
        )

        results_single = simulate(pop_single, params_single)
        final_size_single = length(results_single[end])

        # With the fix, single strain and multiple strain simulations should have
        # similar mortality rates (within tolerance)
        @test abs(final_pop_size - final_size_single) <= tolerance * 0.5  # Smaller tolerance for comparison
    end

    @testset "Strain Introduction Tests" begin
        # Test :simultaneous introduction
        Random.seed!(111)
        pop = Population([Individual(2, 20) for _ in 1:100])

        params_simultaneous = SimulationParameters(
            [SIModel(0.1, 0.0), SIModel(0.1, 0.0)],
            [1.0 1.0; 1.0 1.0],
            0.0, 0.0, 18,
            :simultaneous,  # Both strains introduced simultaneously
            20
        )

        results_sim = simulate(pop, params_simultaneous)

        # Both strains should be present early in simulation
        early_pop = results_sim[2]  # Second time step
        strain1_early = count(ind -> ind[1, 3], early_pop)
        strain2_early = count(ind -> ind[2, 3], early_pop)

        @test strain1_early > 0  # Strain 1 should be introduced
        @test strain2_early > 0  # Strain 2 should be introduced

        # Test :random introduction
        Random.seed!(222)
        pop_random = Population([Individual(3, 20) for _ in 1:100])

        params_random = SimulationParameters(
            [SIModel(0.05, 0.0), SIModel(0.05, 0.0), SIModel(0.05, 0.0)],
            ones(3, 3),
            0.0, 0.0, 18,
            :random,  # Random strain introductions
            50
        )

        results_random = simulate(pop_random, params_random)

        # At least some strains should be introduced during simulation
        final_pop = results_random[end]
        total_infections = sum(count(ind -> ind[i, 3], final_pop) for i in 1:3)
        @test total_infections > 0  # Some infections should occur

        # Test that different strains are introduced at different times (probabilistically)
        strain_first_appearance = zeros(Int, 3)
        for t in 1:length(results_random)
            pop_t = results_random[t]
            for strain in 1:3
                if strain_first_appearance[strain] == 0 && count(ind -> ind[strain, 3], pop_t) > 0
                    strain_first_appearance[strain] = t
                end
            end
        end

        # At least one strain should appear after timestep 1
        @test any(strain_first_appearance .> 1)

        # Test :none introduction (baseline)
        Random.seed!(333)
        pop_none = Population([Individual(2, 20) for _ in 1:100])

        params_none = SimulationParameters(
            [SIModel(0.1, 0.0), SIModel(0.1, 0.0)],
            [1.0 1.0; 1.0 1.0],
            0.0, 0.0, 18,
            :none,
            20
        )

        results_none = simulate(pop_none, params_none)

        # No introductions should occur - only initial infections remain
        all_pops = results_none
        for pop_t in all_pops
            strain1_count = count(ind -> ind[1, 3], pop_t)
            strain2_count = count(ind -> ind[2, 3], pop_t)
            @test strain1_count == 0  # No strain 1 infections (none initially seeded)
            @test strain2_count == 0  # No strain 2 infections (none initially seeded)
        end
    end

    @testset "Aging Process Tests" begin
        # Test that individuals age over time
        Random.seed!(444)
        initial_ages = [5, 10, 15, 20, 25]
        pop = Population([Individual(1, age) for age in initial_ages])

        params = SimulationParameters(
            [SIModel(0.0, 0.0)],  # No disease dynamics
            [1.0;;],
            0.0,  # No mortality
            0.0,  # No births
            100,  # High maturity age
            :none,
            10    # 10 time steps
        )

        results = simulate(pop, params)
        final_pop = results[end]

        # Check that all individuals aged by 9 time steps (10 steps - 1 initial)
        final_ages = [ind.age for ind in final_pop]
        expected_ages = initial_ages .+ 9

        @test length(final_ages) == length(expected_ages)  # No deaths/births
        @test sort(final_ages) == sort(expected_ages)      # Ages increased correctly

        # Test age-dependent reproduction
        Random.seed!(555)
        young_pop = Population([Individual(1, 5) for _ in 1:20])   # Below maturity
        old_pop = Population([Individual(1, 25) for _ in 1:20])    # Above maturity

        params_young = SimulationParameters(
            [SIModel(0.0, 0.0)], [1.0;;], 0.0, 0.1, 18, :none, 5
        )

        params_old = SimulationParameters(
            [SIModel(0.0, 0.0)], [1.0;;], 0.0, 0.1, 18, :none, 5
        )

        results_young = simulate(young_pop, params_young)
        results_old = simulate(old_pop, params_old)

        # Old population should have more births than young population
        final_young_size = length(results_young[end])
        final_old_size = length(results_old[end])

        @test final_old_size > final_young_size  # More reproduction in mature population
        @test final_young_size >= 20            # Young pop should not reproduce much
    end

    @testset "Parameter Validation Tests" begin
        # Test SimulationParameters validation
        valid_models = [SIModel(0.1, 0.01)]
        valid_interactions = [1.0;;]

        # Test invalid number of models vs interaction matrix size
        @test_throws ArgumentError SimulationParameters(
            [SIModel(0.1, 0.01), SIModel(0.1, 0.01)],  # 2 models
            [1.0;;],                                    # 1x1 matrix
            0.01, 0.01, 18, :none, 10
        )

        # Test negative parameters
        @test_throws ArgumentError SimulationParameters(
            valid_models, valid_interactions,
            -0.1, 0.01, 18, :none, 10  # Negative base mortality
        )

        @test_throws ArgumentError SimulationParameters(
            valid_models, valid_interactions,
            0.01, -0.1, 18, :none, 10  # Negative fecundity
        )

        @test_throws ArgumentError SimulationParameters(
            valid_models, valid_interactions,
            0.01, 0.01, -1, :none, 10  # Negative age maturity
        )

        @test_throws ArgumentError SimulationParameters(
            valid_models, valid_interactions,
            0.01, 0.01, 18, :none, 0   # Zero time steps
        )

        # Test invalid introduction type
        @test_throws ArgumentError SimulationParameters(
            valid_models, valid_interactions,
            0.01, 0.01, 18, :invalid, 10
        )

        # Test valid parameters don't throw
        @test_nowarn SimulationParameters(
            valid_models, valid_interactions,
            0.01, 0.01, 18, :none, 10
        )

        @test_nowarn SimulationParameters(
            valid_models, valid_interactions,
            0.01, 0.01, 18, :simultaneous, 10
        )

        @test_nowarn SimulationParameters(
            valid_models, valid_interactions,
            0.01, 0.01, 18, :random, 10
        )
    end

    @testset "Edge Cases and Robustness Tests" begin
        # Test simulation with empty population
        empty_pop = Population(Individual[])
        params = SimulationParameters(
            [SIModel(0.1, 0.01)], [1.0;;], 0.01, 0.01, 18, :none, 5
        )

        results_empty = simulate(empty_pop, params)
        @test length(results_empty) == 5    # Should return 5 time steps
        @test all(length.(results_empty) .== 0)  # All populations should be empty

        # Test simulation with single individual
        single_pop = Population([Individual(1, 20)])
        results_single = simulate(single_pop, params)

        @test length(results_single) == 5
        # Population size could vary due to mortality/reproduction, but should exist
        @test all(length.(results_single) .>= 0)

        # Test with very high mortality (population extinction)
        Random.seed!(666)
        high_mortality_pop = Population([Individual(1, 20) for _ in 1:10])
        high_mortality_params = SimulationParameters(
            [SIModel(0.0, 0.9)],  # Very high disease mortality
            [1.0;;], 0.5, 0.0, 18, :none, 20  # High base mortality too, no reproduction
        )

        # Seed one infection
        high_mortality_pop[1][1, 1] = false
        high_mortality_pop[1][1, 3] = true

        results_extinction = simulate(high_mortality_pop, high_mortality_params)

        # Population should decline significantly or go extinct
        final_size = length(results_extinction[end])
        @test final_size < 10  # Should have fewer individuals than started with

        # Test with very high reproduction
        Random.seed!(777)
        high_repro_pop = Population([Individual(1, 25) for _ in 1:5])  # All mature
        high_repro_params = SimulationParameters(
            [SIModel(0.0, 0.0)], [1.0;;], 0.0, 0.8, 18, :none, 10  # Very high fecundity
        )

        results_growth = simulate(high_repro_pop, high_repro_params)

        # Population should grow substantially
        final_growth_size = length(results_growth[end])
        @test final_growth_size > 10  # Should be much larger than initial

        # Test interaction matrix edge cases
        Random.seed!(888)
        edge_pop = Population([Individual(2, 20) for _ in 1:20])

        # Seed infections
        edge_pop[1][1, 1] = false; edge_pop[1][1, 3] = true
        edge_pop[2][2, 1] = false; edge_pop[2][2, 3] = true

        # Test with extreme facilitation
        extreme_facilitation = [1.0 5.0; 5.0 1.0]  # Very strong facilitation
        params_facilitation = SimulationParameters(
            [SIModel(0.1, 0.0), SIModel(0.1, 0.0)],
            extreme_facilitation, 0.0, 0.0, 18, :none, 15
        )

        results_facilitation = simulate(edge_pop, params_facilitation)
        @test length(results_facilitation) == 15  # Should complete without error

        # Test with extreme competition
        Random.seed!(888)  # Same seed for comparison
        edge_pop2 = Population([Individual(2, 20) for _ in 1:20])
        edge_pop2[1][1, 1] = false; edge_pop2[1][1, 3] = true
        edge_pop2[2][2, 1] = false; edge_pop2[2][2, 3] = true

        extreme_competition = [1.0 0.1; 0.1 1.0]   # Very strong competition
        params_competition = SimulationParameters(
            [SIModel(0.1, 0.0), SIModel(0.1, 0.0)],
            extreme_competition, 0.0, 0.0, 18, :none, 15
        )

        results_competition = simulate(edge_pop2, params_competition)
        @test length(results_competition) == 15  # Should complete without error

        # With extreme competition, final infections should be lower than with facilitation
        final_facilitation_infections = sum(count(ind -> ind[i, 3], results_facilitation[end]) for i in 1:2)
        final_competition_infections = sum(count(ind -> ind[i, 3], results_competition[end]) for i in 1:2)

        # This is probabilistic, but competition should generally reduce spread
        @test final_competition_infections <= final_facilitation_infections
    end

    @testset "Disease State Transitions" begin
        # Test that disease state transitions work correctly for each model type

        # Test SEIR latency period
        Random.seed!(999)
        seir_pop = Population([Individual(1, 20) for _ in 1:100])

        # Expose some individuals initially (not infected)
        for i in 1:20
            seir_pop[i][1, 1] = false  # Not susceptible
            seir_pop[i][1, 2] = true   # Exposed
        end

        latency_period = 5
        seir_params = SimulationParameters(
            [SEIRModel(0.0, 0.0, 0.0, latency_period)],  # No transmission/mortality/recovery
            [1.0;;], 0.0, 0.0, 100, :none, latency_period * 2
        )

        results_seir = simulate(seir_pop, seir_params)

        # Check that exposed individuals transition to infected over time
        exposed_counts = [count(ind -> ind[1, 2], pop) for pop in results_seir]
        infected_counts = [count(ind -> ind[1, 3], pop) for pop in results_seir]

        # Initially should have exposed individuals
        @test exposed_counts[1] == 20
        @test infected_counts[1] == 0

        # Over time, some exposed should become infected
        @test infected_counts[end] > 0
        @test exposed_counts[end] < 20

        # Test SEIRS immunity loss
        Random.seed!(42)
        seirs_pop = Population([Individual(1, 20) for _ in 1:100])

        # Make some individuals recovered initially
        for i in 1:30
            seirs_pop[i][1, 1] = false  # Not susceptible
            seirs_pop[i][1, 4] = true   # Recovered
        end

        immunity_loss_rate = 0.3
        seirs_params = SimulationParameters(
            [SEIRSModel(0.0, 0.0, 0.0, 3, immunity_loss_rate)],
            [1.0;;], 0.0, 0.0, 100, :none, 20
        )

        results_seirs = simulate(seirs_pop, seirs_params)

        # Check immunity loss over time
        recovered_counts = [count(ind -> ind[1, 4], pop) for pop in results_seirs]
        susceptible_counts = [count(ind -> ind[1, 1], pop) for pop in results_seirs]

        # Initially should have 30 recovered, rest susceptible
        @test recovered_counts[1] == 30
        @test susceptible_counts[1] == 70

        # Over time, some recovered should lose immunity and become susceptible again
        @test susceptible_counts[end] > 70  # Some recovered individuals lost immunity
    end

    @testset "Population Methods and Indexing" begin
        # Test Population indexing and methods
        ind1 = Individual(2, 10)
        ind2 = Individual(2, 20)
        ind3 = Individual(2, 30)

        pop = Population([ind1, ind2, ind3])

        # Test indexing
        @test pop[1] === ind1
        @test pop[2] === ind2
        @test pop[3] === ind3

        # Test length
        @test length(pop) == 3

        # Test iteration
        collected = collect(pop)
        @test length(collected) == 3
        @test collected[1] === ind1

        # Test n_strains helper
        @test CoinfectionSimulator.n_strains(pop) == 2

        # Test add_individual! function
        ind4 = Individual(2, 40)
        CoinfectionSimulator.add_individual!(pop, ind4)
        @test length(pop) == 4
        @test pop[4] === ind4

        # Test error when adding individual with different strain count
        ind_wrong_strains = Individual(3, 50)
        @test_throws ArgumentError CoinfectionSimulator.add_individual!(pop, ind_wrong_strains)

        # Test empty population
        empty_pop = Population(Individual[])
        @test length(empty_pop) == 0

        # First individual sets strain count for empty population
        CoinfectionSimulator.add_individual!(empty_pop, Individual(4, 25))
        @test length(empty_pop) == 1
        @test CoinfectionSimulator.n_strains(empty_pop) == 4
    end

    @testset "Simulation Helper Functions" begin
        # Test introduce_infections! function
        Random.seed!(1111)
        test_pop = Population([Individual(3, 20) for _ in 1:50])

        # All start susceptible - simulate introducing strain 2 at timestep 5
        intro_schedule = [0, 0, 5]  # Strain 3 introduced at timestep 5

        # Before introduction
        initial_strain3_count = count(ind -> ind[3, 3], test_pop)
        @test initial_strain3_count == 0

        # Introduce infections
        CoinfectionSimulator.introduce_infections!(test_pop, 5, intro_schedule)

        # After introduction, should have exactly 1 infection of strain 3
        final_strain3_count = count(ind -> ind[3, 3], test_pop)
        @test final_strain3_count == 1

        # Test introducing multiple strains simultaneously
        multi_pop = Population([Individual(4, 20) for _ in 1:100])
        multi_schedule = [1, 1, 0, 1]  # Strains 1, 2, and 4 at timestep 1

        CoinfectionSimulator.introduce_infections!(multi_pop, 1, multi_schedule)

        strain1_count = count(ind -> ind[1, 3], multi_pop)
        strain2_count = count(ind -> ind[2, 3], multi_pop)
        strain3_count = count(ind -> ind[3, 3], multi_pop)
        strain4_count = count(ind -> ind[4, 3], multi_pop)

        @test strain1_count == 1
        @test strain2_count == 1
        @test strain3_count == 0  # Not introduced
        @test strain4_count == 1

        # Test no introduction when timestep doesn't match
        no_intro_pop = Population([Individual(2, 20) for _ in 1:20])
        CoinfectionSimulator.introduce_infections!(no_intro_pop, 3, [5, 5])  # Wrong timestep

        @test count(ind -> ind[1, 3], no_intro_pop) == 0
        @test count(ind -> ind[2, 3], no_intro_pop) == 0
    end

    @testset "Breeding and Demographics" begin
        # Test breeding! function directly
        Random.seed!(1212)
        breeding_pop = Population([Individual(2, 25) for _ in 1:20])  # All mature

        initial_size = length(breeding_pop)
        CoinfectionSimulator.breeding!(breeding_pop, 0.3, 18)  # 30% fecundity, mature at 18

        @test length(breeding_pop) >= initial_size  # Should have births

        # Test that newborns are age 0 and all susceptible
        ages = [ind.age for ind in breeding_pop]
        newborn_count = count(age -> age == 0, ages)

        if newborn_count > 0
            newborns = filter(ind -> ind.age == 0, breeding_pop.individuals)
            for newborn in newborns
                @test all(newborn.state[:, 1])  # All strains susceptible
                @test !any(newborn.state[:, 2:4])  # Not in other states
            end
        end

        # Test no breeding with zero fecundity
        no_breed_pop = Population([Individual(1, 30) for _ in 1:10])
        initial_no_breed_size = length(no_breed_pop)
        CoinfectionSimulator.breeding!(no_breed_pop, 0.0, 18)  # Zero fecundity

        @test length(no_breed_pop) == initial_no_breed_size  # No change

        # Test no breeding with immature population
        young_pop = Population([Individual(1, 10) for _ in 1:20])  # All below maturity
        initial_young_size = length(young_pop)
        CoinfectionSimulator.breeding!(young_pop, 0.5, 18)  # High fecundity but immature

        @test length(young_pop) == initial_young_size  # No breeding
    end

    @testset "Disease Process Functions" begin
        # Test individual disease process functions
        Random.seed!(1313)

        # Test process_recovery!
        recovery_pop = Population([Individual(1, 20) for _ in 1:100])

        # Make half infected
        for i in 1:50
            recovery_pop[i][1, 1] = false  # Not susceptible
            recovery_pop[i][1, 3] = true   # Infected
        end

        initial_infected = count(ind -> ind[1, 3], recovery_pop)
        @test initial_infected == 50

        # Apply recovery with high rate
        alive = fill(true, length(recovery_pop))
        CoinfectionSimulator.process_recovery!(recovery_pop, alive, 1, 0.8)  # 80% recovery rate

        final_infected = count(ind -> ind[1, 3], recovery_pop)
        final_recovered = count(ind -> ind[1, 4], recovery_pop)

        @test final_infected < 50  # Some should recover
        @test final_recovered > 0  # Some should be recovered
        @test final_infected + final_recovered == 50  # Conservation

        # Test process_immunity_loss!
        immunity_pop = Population([Individual(1, 20) for _ in 1:100])

        # Make all recovered
        for i in 1:100
            immunity_pop[i][1, 1] = false  # Not susceptible
            immunity_pop[i][1, 4] = true   # Recovered
        end

        initial_recovered = count(ind -> ind[1, 4], immunity_pop)
        @test initial_recovered == 100

        alive = fill(true, length(immunity_pop))
        CoinfectionSimulator.process_immunity_loss!(immunity_pop, alive, 1, 0.3)  # 30% immunity loss

        final_recovered = count(ind -> ind[1, 4], immunity_pop)
        final_susceptible = count(ind -> ind[1, 1], immunity_pop)

        @test final_recovered < 100  # Some should lose immunity
        @test final_susceptible > 0  # Some should become susceptible again
        @test final_recovered + final_susceptible == 100  # Conservation

        # Test process_latent_infections! (exposed to infected)
        latent_pop = Population([Individual(1, 20) for _ in 1:100])

        # Make all exposed
        for i in 1:100
            latent_pop[i][1, 1] = false  # Not susceptible
            latent_pop[i][1, 2] = true   # Exposed
        end

        initial_exposed = count(ind -> ind[1, 2], latent_pop)
        @test initial_exposed == 100

        alive = fill(true, length(latent_pop))
        CoinfectionSimulator.process_latent_infections!(latent_pop, alive, 1, 2)  # Latency period = 2

        final_exposed = count(ind -> ind[1, 2], latent_pop)
        final_infected = count(ind -> ind[1, 3], latent_pop)

        @test final_exposed < 100  # Some should become infected
        @test final_infected > 0   # Some should be infected
        @test final_exposed + final_infected == 100  # Conservation
    end

    @testset "Mortality Processes" begin
        # Test apply_base_mortality!
        Random.seed!(1414)
        mortality_pop = Population([Individual(1, 20) for _ in 1:1000])  # Large population for statistics

        alive = fill(true, length(mortality_pop))
        mortality_set = Set{Int}()

        # Apply moderate base mortality
        CoinfectionSimulator.apply_base_mortality!(mortality_pop, alive, 0.05, mortality_set)

        # Should have some deaths
        @test length(mortality_set) > 0
        @test length(mortality_set) < 1000  # But not everyone

        # Test disease mortality
        disease_mortality_pop = Population([Individual(1, 20) for _ in 1:100])

        # Make all infected
        for i in 1:100
            disease_mortality_pop[i][1, 1] = false
            disease_mortality_pop[i][1, 3] = true
        end

        alive = fill(true, length(disease_mortality_pop))
        mortality_set = Set{Int}()

        CoinfectionSimulator.process_disease_mortality!(disease_mortality_pop, alive, 1, 0.1, mortality_set)

        # Should have some disease-related deaths
        @test length(mortality_set) > 0

        # All deaths should be from infected individuals
        for dead_idx in mortality_set
            @test disease_mortality_pop[dead_idx][1, 3]  # Should be infected
        end
    end

    @testset "Comprehensive Integration Tests" begin
        # Test full simulation with all components working together
        Random.seed!(1515)

        # Create a complex scenario with multiple strains, demographics, and introductions
        complex_pop = Population([Individual(4, rand(18:50)) for _ in 1:200])

        # Mix of disease models
        complex_models = [
            SIModel(0.1, 0.01),      # Strain 1: SI
            SIRModel(0.12, 0.01, 0.15),  # Strain 2: SIR
            SEIRModel(0.08, 0.005, 0.12, 4),  # Strain 3: SEIR
            SEIRSModel(0.09, 0.008, 0.1, 3, 0.05)  # Strain 4: SEIRS
        ]

        # Complex interaction matrix
        complex_interactions = [
            1.0  0.8  1.2  0.9;
            1.1  1.0  0.7  1.3;
            0.9  1.2  1.0  1.1;
            1.4  0.8  1.1  1.0
        ]

        complex_params = SimulationParameters(
            complex_models,
            complex_interactions,
            0.005,   # Base mortality
            0.02,    # Fecundity
            18,      # Age maturity
            :random, # Random strain introductions
            100      # Long simulation
        )

        # Should complete without error
        @test_nowarn complex_results = simulate(complex_pop, complex_params)

        complex_results = simulate(complex_pop, complex_params)

        # Basic sanity checks
        @test length(complex_results) == 100  # 100 time steps
        @test all(length.(complex_results) .>= 0)  # Valid population sizes

        # Should have some disease activity
        final_pop = complex_results[end]
        if length(final_pop) > 0
            total_infections = sum(count(ind -> ind[strain, 3], final_pop) for strain in 1:4)
            total_recovered = sum(count(ind -> ind[strain, 4], final_pop) for strain in 1:4)

            # Should have had some disease activity over 100 timesteps
            @test (total_infections + total_recovered) > 0
        end

        # Ages should have increased
        if length(final_pop) > 0
            ages = [ind.age for ind in final_pop]
            # At least some individuals should be older than initial range
            @test any(age -> age > 50, ages) || any(age -> age < 18, ages)  # Either aging or new births
        end
    end

end
