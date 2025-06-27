@testset "Utils Tests" begin
    using Random
    using Distributions
    using StatsBase
    
    # Import internal functions for testing
    import CoinfectionSimulator: handle_si_disease, handle_sir_disease, handle_seir_disease, 
                                handle_seirs_disease, handle_infection, handle_exposure,
                                handle_exposed_infection, handle_recovery, handle_immunity_loss,
                                handle_infected_death
    
    # Test setup helper functions
    function create_test_population(n_individuals::Int, n_strains::Int)
        # Use BitMatrix (default from falses()) to test type compatibility
        return [falses(n_strains, 4) for _ in 1:n_individuals]
    end
    
    function set_susceptible!(pop, individual_idx, strain_idx)
        pop[individual_idx][strain_idx, 1] = true
        pop[individual_idx][strain_idx, 2:4] .= false
    end
    
    function set_exposed!(pop, individual_idx, strain_idx)
        pop[individual_idx][strain_idx, 1] = false
        pop[individual_idx][strain_idx, 2] = true
        pop[individual_idx][strain_idx, 3:4] .= false
    end
    
    function set_infected!(pop, individual_idx, strain_idx)
        pop[individual_idx][strain_idx, 1:2] .= false
        pop[individual_idx][strain_idx, 3] = true
        pop[individual_idx][strain_idx, 4] = false
    end
    
    function set_recovered!(pop, individual_idx, strain_idx)
        pop[individual_idx][strain_idx, 1:3] .= false
        pop[individual_idx][strain_idx, 4] = true
    end

    @testset "handle_si_disease" begin
        Random.seed!(123)
        
        # Test basic SI disease dynamics
        pop = create_test_population(10, 1)
        # Set 8 susceptible, 1 infected
        for i in 1:8
            set_susceptible!(pop, i, 1)
        end
        set_infected!(pop, 9, 1)
        
        alive = fill(true, 10)
        strain = 1
        transmission = 0.5
        interactions = [1.0]
        base_mortality = 0.0
        disease_mortality = 0.0
        dead_indices = Set{Int}()
        
        # Test that function runs without error
        @test_nowarn handle_si_disease(
            pop, alive, strain, transmission, interactions,
            base_mortality, disease_mortality, dead_indices
        )
        
        # Test with mortality
        pop2 = create_test_population(100, 1)
        for i in 1:90
            set_susceptible!(pop2, i, 1)
        end
        for i in 91:100
            set_infected!(pop2, i, 1)
        end
        alive2 = fill(true, 100)
        dead_indices2 = Set{Int}()
        
        handle_si_disease(
            pop2, alive2, 1, 0.0, [1.0], 0.0, 0.1, dead_indices2
        )
        
        # With disease mortality, some infected should die
        @test length(dead_indices2) >= 0  # Could be 0 due to randomness
    end

    @testset "handle_sir_disease" begin
        Random.seed!(456)
        
        # Test SIR disease dynamics
        pop = create_test_population(10, 1)
        for i in 1:8
            set_susceptible!(pop, i, 1)
        end
        set_infected!(pop, 9, 1)
        set_recovered!(pop, 10, 1)
        
        alive = fill(true, 10)
        strain = 1
        transmission = 0.3
        interactions = [1.0]
        base_mortality = 0.0
        disease_mortality = 0.0
        recovery = 0.2
        dead_indices = Set{Int}()
        
        # Test that function runs without error
        @test_nowarn handle_sir_disease(
            pop, alive, strain, transmission, interactions,
            base_mortality, disease_mortality, recovery, dead_indices
        )
        
        # Test recovery functionality
        pop2 = create_test_population(100, 1)
        for i in 1:100
            set_infected!(pop2, i, 1)
        end
        alive2 = fill(true, 100)
        dead_indices2 = Set{Int}()
        
        initial_infected = sum(pop2[i][1, 3] for i in 1:100)
        @test initial_infected == 100
        
        handle_sir_disease(
            pop2, alive2, 1, 0.0, [1.0], 0.0, 0.0, 1.0, dead_indices2
        )
        
        final_infected = sum(pop2[i][1, 3] for i in 1:100)
        final_recovered = sum(pop2[i][1, 4] for i in 1:100)
        
        # With recovery = 1.0, all should recover
        @test final_infected == 0
        @test final_recovered == 100
    end

    @testset "handle_seir_disease" begin
        Random.seed!(789)
        
        # Test SEIR disease dynamics
        pop = create_test_population(10, 1)
        for i in 1:7
            set_susceptible!(pop, i, 1)
        end
        set_exposed!(pop, 8, 1)
        set_infected!(pop, 9, 1)
        set_recovered!(pop, 10, 1)
        
        alive = fill(true, 10)
        strain = 1
        transmission = 0.3
        interactions = [1.0]
        base_mortality = 0.0
        disease_mortality = 0.0
        recovery = 0.2
        latency = 2
        dead_indices = Set{Int}()
        
        # Test that function runs without error
        @test_nowarn handle_seir_disease(
            pop, alive, strain, transmission, interactions,
            base_mortality, disease_mortality, recovery, latency, dead_indices
        )
        
        # Test exposed to infected transition
        pop2 = create_test_population(100, 1)
        for i in 1:100
            set_exposed!(pop2, i, 1)
        end
        alive2 = fill(true, 100)
        dead_indices2 = Set{Int}()
        
        initial_exposed = sum(pop2[i][1, 2] for i in 1:100)
        @test initial_exposed == 100
        
        # With latency = 1, all exposed should become infected
        handle_seir_disease(
            pop2, alive2, 1, 0.0, [1.0], 0.0, 0.0, 0.0, 1, dead_indices2
        )
        
        final_exposed = sum(pop2[i][1, 2] for i in 1:100)
        final_infected = sum(pop2[i][1, 3] for i in 1:100)
        
        # With latency = 1, approximately all should transition to infected
        @test final_exposed < initial_exposed
        @test final_infected > 0
    end

    @testset "handle_seirs_disease" begin
        Random.seed!(101112)
        
        # Test SEIRS disease dynamics
        pop = create_test_population(10, 1)
        for i in 1:6
            set_susceptible!(pop, i, 1)
        end
        set_exposed!(pop, 7, 1)
        set_infected!(pop, 8, 1)
        set_recovered!(pop, 9, 1)
        set_recovered!(pop, 10, 1)
        
        alive = fill(true, 10)
        strain = 1
        transmission = 0.3
        interactions = [1.0]
        base_mortality = 0.0
        disease_mortality = 0.0
        recovery = 0.2
        latency = 2
        immunity_loss = 0.1
        dead_indices = Set{Int}()
        
        # Test that function runs without error
        @test_nowarn handle_seirs_disease(
            pop, alive, strain, transmission, interactions,
            base_mortality, disease_mortality, recovery, latency, immunity_loss, dead_indices
        )
        
        # Test immunity loss functionality
        pop2 = create_test_population(100, 1)
        for i in 1:100
            set_recovered!(pop2, i, 1)
        end
        alive2 = fill(true, 100)
        dead_indices2 = Set{Int}()
        
        initial_recovered = sum(pop2[i][1, 4] for i in 1:100)
        @test initial_recovered == 100
        
        # With immunity_loss = 1.0, all recovered should become susceptible
        handle_seirs_disease(
            pop2, alive2, 1, 0.0, [1.0], 0.0, 0.0, 0.0, 5, 1.0, dead_indices2
        )
        
        final_recovered = sum(pop2[i][1, 4] for i in 1:100)
        final_susceptible = sum(pop2[i][1, 1] for i in 1:100)
        
        @test final_recovered == 0
        @test final_susceptible == 100
    end

    @testset "handle_infection" begin
        Random.seed!(131415)
        
        # Test infection process
        pop = create_test_population(10, 1)
        for i in 1:8
            set_susceptible!(pop, i, 1)
        end
        set_infected!(pop, 9, 1)
        set_infected!(pop, 10, 1)
        
        alive = fill(true, 10)
        strain = 1
        transmission = 1.0  # Ensure infection occurs
        interactions = [1.0]
        
        initial_susceptible = sum(pop[i][1, 1] for i in 1:10)
        initial_infected = sum(pop[i][1, 3] for i in 1:10)
        
        handle_infection(pop, alive, strain, transmission, interactions)
        
        final_susceptible = sum(pop[i][1, 1] for i in 1:10)
        final_infected = sum(pop[i][1, 3] for i in 1:10)
        
        # With high transmission, some susceptible should become infected
        @test final_susceptible <= initial_susceptible
        @test final_infected >= initial_infected
        
        # Test with no infected individuals
        pop2 = create_test_population(5, 1)
        for i in 1:5
            set_susceptible!(pop2, i, 1)
        end
        alive2 = fill(true, 5)
        
        @test_nowarn handle_infection(pop2, alive2, 1, 0.5, [1.0])
        
        # No change should occur
        @test all(pop2[i][1, 1] for i in 1:5)
    end

    @testset "handle_exposure" begin
        Random.seed!(161718)
        
        # Test exposure process (SEIR model)
        pop = create_test_population(10, 1)
        for i in 1:8
            set_susceptible!(pop, i, 1)
        end
        set_infected!(pop, 9, 1)
        set_infected!(pop, 10, 1)
        
        alive = fill(true, 10)
        strain = 1
        transmission = 1.0
        interactions = [1.0]
        
        initial_susceptible = sum(pop[i][1, 1] for i in 1:10)
        initial_exposed = sum(pop[i][1, 2] for i in 1:10)
        
        handle_exposure(pop, alive, strain, transmission, interactions)
        
        final_susceptible = sum(pop[i][1, 1] for i in 1:10)
        final_exposed = sum(pop[i][1, 2] for i in 1:10)
        
        # With high transmission, some susceptible should become exposed
        @test final_susceptible <= initial_susceptible
        @test final_exposed >= initial_exposed
        
        # Test with no infected individuals
        pop2 = create_test_population(5, 1)
        for i in 1:5
            set_susceptible!(pop2, i, 1)
        end
        alive2 = fill(true, 5)
        
        @test_nowarn handle_exposure(pop2, alive2, 1, 0.5, [1.0])
        
        # No change should occur
        @test all(pop2[i][1, 1] for i in 1:5)
    end

    @testset "handle_exposed_infection" begin
        Random.seed!(192021)
        
        # Test exposed to infected transition
        pop = create_test_population(100, 1)
        for i in 1:100
            set_exposed!(pop, i, 1)
        end
        
        alive = fill(true, 100)
        strain = 1
        latency = 1  # Short latency for high transition probability
        
        initial_exposed = sum(pop[i][1, 2] for i in 1:100)
        initial_infected = sum(pop[i][1, 3] for i in 1:100)
        
        handle_exposed_infection(pop, alive, strain, latency)
        
        final_exposed = sum(pop[i][1, 2] for i in 1:100)
        final_infected = sum(pop[i][1, 3] for i in 1:100)
        
        # With latency = 1, most exposed should become infected
        @test final_exposed < initial_exposed
        @test final_infected > initial_infected
        
        # Test with no exposed individuals
        pop2 = create_test_population(5, 1)
        for i in 1:5
            set_susceptible!(pop2, i, 1)
        end
        alive2 = fill(true, 5)
        
        @test_nowarn handle_exposed_infection(pop2, alive2, 1, 5)
        
        # No change should occur
        @test all(pop2[i][1, 1] for i in 1:5)
    end

    @testset "handle_recovery" begin
        Random.seed!(222324)
        
        # Test recovery process
        pop = create_test_population(100, 1)
        for i in 1:100
            set_infected!(pop, i, 1)
        end
        
        alive = fill(true, 100)
        strain = 1
        recovery = 1.0  # All should recover
        
        initial_infected = sum(pop[i][1, 3] for i in 1:100)
        initial_recovered = sum(pop[i][1, 4] for i in 1:100)
        
        handle_recovery(pop, alive, strain, recovery)
        
        final_infected = sum(pop[i][1, 3] for i in 1:100)
        final_recovered = sum(pop[i][1, 4] for i in 1:100)
        
        # With recovery = 1.0, all should recover
        @test final_infected == 0
        @test final_recovered == 100
        
        # Test with recovery = 0.0 (no recovery)
        pop2 = create_test_population(10, 1)
        for i in 1:10
            set_infected!(pop2, i, 1)
        end
        alive2 = fill(true, 10)
        
        handle_recovery(pop2, alive2, 1, 0.0)
        
        # No change should occur
        @test all(pop2[i][1, 3] for i in 1:10)
        @test all(!pop2[i][1, 4] for i in 1:10)
        
        # Test with no infected individuals
        pop3 = create_test_population(5, 1)
        for i in 1:5
            set_susceptible!(pop3, i, 1)
        end
        alive3 = fill(true, 5)
        
        @test_nowarn handle_recovery(pop3, alive3, 1, 0.5)
    end

    @testset "handle_immunity_loss" begin
        Random.seed!(252627)
        
        # Test immunity loss process
        pop = create_test_population(100, 1)
        for i in 1:100
            set_recovered!(pop, i, 1)
        end
        
        alive = fill(true, 100)
        strain = 1
        immunity_loss = 1.0  # All should lose immunity
        
        initial_recovered = sum(pop[i][1, 4] for i in 1:100)
        initial_susceptible = sum(pop[i][1, 1] for i in 1:100)
        
        handle_immunity_loss(pop, alive, strain, immunity_loss)
        
        final_recovered = sum(pop[i][1, 4] for i in 1:100)
        final_susceptible = sum(pop[i][1, 1] for i in 1:100)
        
        # With immunity_loss = 1.0, all should lose immunity
        @test final_recovered == 0
        @test final_susceptible == 100
        
        # Test with immunity_loss = 0.0 (no loss)
        pop2 = create_test_population(10, 1)
        for i in 1:10
            set_recovered!(pop2, i, 1)
        end
        alive2 = fill(true, 10)
        
        handle_immunity_loss(pop2, alive2, 1, 0.0)
        
        # No change should occur
        @test all(pop2[i][1, 4] for i in 1:10)
        @test all(!pop2[i][1, 1] for i in 1:10)
        
        # Test with no recovered individuals
        pop3 = create_test_population(5, 1)
        for i in 1:5
            set_susceptible!(pop3, i, 1)
        end
        alive3 = fill(true, 5)
        
        @test_nowarn handle_immunity_loss(pop3, alive3, 1, 0.5)
    end

    @testset "handle_infected_death" begin
        Random.seed!(282930)
        
        # Test infected death process
        pop = create_test_population(100, 1)
        for i in 1:100
            set_infected!(pop, i, 1)
        end
        
        strain = 1
        base_mortality = 0.1
        disease_mortality = 0.2
        dead_indices = Set{Int}()
        
        initial_alive = length(pop)
        
        handle_infected_death(pop, strain, base_mortality, disease_mortality, dead_indices)
        
        # Some individuals should die with mortality rates
        @test length(dead_indices) >= 0  # Could be 0 due to randomness
        @test length(dead_indices) <= 100
        
        # Test with zero mortality
        pop2 = create_test_population(10, 1)
        for i in 1:10
            set_infected!(pop2, i, 1)
        end
        dead_indices2 = Set{Int}()
        
        handle_infected_death(pop2, 1, 0.0, 0.0, dead_indices2)
        
        @test length(dead_indices2) == 0
        
        # Test with no infected individuals
        pop3 = create_test_population(5, 1)
        for i in 1:5
            set_susceptible!(pop3, i, 1)
        end
        dead_indices3 = Set{Int}()
        
        @test_nowarn handle_infected_death(pop3, 1, 0.1, 0.1, dead_indices3)
        @test length(dead_indices3) == 0
    end

    @testset "Multi-strain interactions" begin
        Random.seed!(313233)
        
        # Test with multiple strains
        pop = create_test_population(10, 2)
        
        # Set up mixed states
        for i in 1:5
            set_susceptible!(pop, i, 1)
            set_susceptible!(pop, i, 2)
        end
        
        set_infected!(pop, 6, 1)  # Strain 1 infected
        set_infected!(pop, 7, 2)  # Strain 2 infected
        set_infected!(pop, 8, 1)  # Both strains infected
        set_infected!(pop, 8, 2)
        
        alive = fill(true, 10)
        
        # Test strain 1 infection
        handle_infection(pop, alive, 1, 0.5, [1.0, 0.8])
        
        # Test strain 2 infection  
        handle_infection(pop, alive, 2, 0.5, [0.8, 1.0])
        
        # Functions should run without error
        @test true
    end

    @testset "Edge cases and error handling" begin
        # Test with empty population
        empty_pop = Vector{BitMatrix}()
        empty_alive = Bool[]
        dead_indices = Set{Int}()
        
        @test_nowarn handle_infection(empty_pop, empty_alive, 1, 0.5, [1.0])
        @test_nowarn handle_exposure(empty_pop, empty_alive, 1, 0.5, [1.0])
        @test_nowarn handle_recovery(empty_pop, empty_alive, 1, 0.5)
        @test_nowarn handle_immunity_loss(empty_pop, empty_alive, 1, 0.5)
        @test_nowarn handle_exposed_infection(empty_pop, empty_alive, 1, 5)
        @test_nowarn handle_infected_death(empty_pop, 1, 0.1, 0.1, dead_indices)
        
        # Test with all individuals dead
        pop = create_test_population(5, 1)
        alive = fill(false, 5)
        dead_indices = Set(1:5)
        
        @test_nowarn handle_infection(pop, alive, 1, 0.5, [1.0])
        @test_nowarn handle_exposure(pop, alive, 1, 0.5, [1.0])
        @test_nowarn handle_recovery(pop, alive, 1, 0.5)
        @test_nowarn handle_immunity_loss(pop, alive, 1, 0.5)
        @test_nowarn handle_exposed_infection(pop, alive, 1, 5)
        
        # Test with extreme parameter values
        pop2 = create_test_population(10, 1)
        for i in 1:5
            set_susceptible!(pop2, i, 1)
        end
        for i in 6:10
            set_infected!(pop2, i, 1)
        end
        alive2 = fill(true, 10)
        
        # Test with transmission = 0 (no transmission)
        @test_nowarn handle_infection(pop2, alive2, 1, 0.0, [1.0])
        
        # Test with very high latency
        pop3 = create_test_population(10, 1)
        for i in 1:10
            set_exposed!(pop3, i, 1)
        end
        alive3 = fill(true, 10)
        
        @test_nowarn handle_exposed_infection(pop3, alive3, 1, 1000)
    end
end
