@testset "Population Sampling Tests" begin
    @testset "Basic Sampling" begin
        # Create a test population with known infection status
        Random.seed!(123)

        # Create a series of populations across 3 time steps
        populations = Vector{Population}(undef, 3)

        # Time step 1: 10 individuals, 2 strains, 2 infected with strain 1, 1 infected with strain 2
        inds1 = Vector{Individual}(undef, 10)
        for i in 1:10
            if i <= 2
                # First two individuals infected with strain 1
                state = BitMatrix([false false true false; true false false false])
            elseif i == 3
                # Third individual infected with strain 2
                state = BitMatrix([true false false false; false false true false])
            else
                # Others all susceptible
                state = BitMatrix([true false false false; true false false false])
            end
            inds1[i] = Individual(state, 10)
        end
        populations[1] = Population(inds1)

        # Time step 2: Same but with more infections
        inds2 = Vector{Individual}(undef, 10)
        for i in 1:10
            if i <= 4
                # Four individuals infected with strain 1
                state = BitMatrix([false false true false; true false false false])
            elseif i <= 6
                # Two individuals infected with strain 2
                state = BitMatrix([true false false false; false false true false])
            elseif i == 7
                # One individual co-infected
                state = BitMatrix([false false true false; false false true false])
            else
                # Others all susceptible
                state = BitMatrix([true false false false; true false false false])
            end
            inds2[i] = Individual(state, 10)
        end
        populations[2] = Population(inds2)

        # Time step 3: Even more infections
        inds3 = Vector{Individual}(undef, 10)
        for i in 1:10
            if i <= 5
                # Five individuals infected with strain 1
                state = BitMatrix([false false true false; true false false false])
            elseif i <= 8
                # Three individuals infected with strain 2
                state = BitMatrix([true false false false; false false true false])
            else
                # Two individuals co-infected
                state = BitMatrix([false false true false; false false true false])
            end
            inds3[i] = Individual(state, 10)
        end
        populations[3] = Population(inds3)

        # Perfect sampling (100% sampled, no errors)
        params_perfect = SamplingParameters(1.0, 0.0, 0.0)
        results_perfect = sample_populations(populations, params_perfect)

        @test size(results_perfect) == (3, 2)  # 3 timesteps, 2 strains

        # Check perfect detection
        @test results_perfect[1, 1] == true  # Strain 1 detected at timestep 1
        @test results_perfect[1, 2] == true  # Strain 2 detected at timestep 1
        @test results_perfect[2, 1] == true  # Strain 1 detected at timestep 2
        @test results_perfect[2, 2] == true  # Strain 2 detected at timestep 2
        @test results_perfect[3, 1] == true  # Strain 1 detected at timestep 3
        @test results_perfect[3, 2] == true  # Strain 2 detected at timestep 3

        # Partial sampling with errors
        # This is somewhat probabilistic, so we use a fixed seed
        Random.seed!(456)
        params_partial = SamplingParameters(0.3, 0.05, 0.2)
        results_partial = sample_populations(populations, params_partial)

        @test size(results_partial) == (3, 2)

        # With partial sampling, results are probabilistic, but we can verify
        # that sampling at least detected strains that were common in the population
        @test any(results_partial[2:3, 1]) == true  # Strain 1 should be detected in at least one later timestep

        # Test with empty population
        empty_pop = Population(Individual[])
        populations_empty = [empty_pop, empty_pop]
        results_empty = sample_populations(populations_empty, params_perfect)
        @test size(results_empty) == (2, 0)  # 2 timesteps, 0 strains
    end

    @testset "Parameter Validation" begin
        # Test invalid parameters
        @test_throws ArgumentError SamplingParameters(-0.1, 0.0, 0.0)  # Negative proportion
        @test_throws ArgumentError SamplingParameters(1.1, 0.0, 0.0)   # Proportion > 1
        @test_throws ArgumentError SamplingParameters(0.5, -0.1, 0.0)  # Negative false positive
        @test_throws ArgumentError SamplingParameters(0.5, 1.1, 0.0)   # False positive > 1
        @test_throws ArgumentError SamplingParameters(0.5, 0.0, -0.1)  # Negative false negative
        @test_throws ArgumentError SamplingParameters(0.5, 0.0, 1.1)   # False negative > 1
    end

    @testset "False Positive and Negative Rates" begin
        # Create a controlled population to test false positive/negative
        Random.seed!(789)

        # Create a population where all individuals are infected with strain 1 only
        all_strain1 = Population([
            Individual(BitMatrix([false false true false; true false false false]), 1)
            for _ in 1:100
        ])

        # Create a population where no individuals are infected
        none_infected = Population([
            Individual(BitMatrix([true false false false; true false false false]), 1)
            for _ in 1:100
        ])

        # Test false negative rate
        params_fn = SamplingParameters(1.0, 0.0, 0.3)  # 30% false negative, no false positive
        result_fn = sample_populations([all_strain1], params_fn)
        @test result_fn[1, 1] == true  # Should still detect strain 1 despite false negatives
        @test result_fn[1, 2] == false # Should not detect strain 2

        # Test false positive rate
        Random.seed!(790)  # Set different seed for reproducibility
        params_fp = SamplingParameters(1.0, 0.3, 0.0)  # 30% false positive, no false negative
        result_fp = sample_populations([none_infected], params_fp)
        # With 100 individuals and 30% false positive rate, we should detect some false positives
        @test result_fp[1, 1] || result_fp[1, 2]
    end


end
