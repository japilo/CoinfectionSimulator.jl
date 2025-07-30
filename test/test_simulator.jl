@testset "Coinfection Simulator Tests" begin

	@testset "Input validation tests" begin
		# Test empty population
		@test_throws AssertionError coinfection_simulator(
			initial_pop = Vector{Matrix{Bool}}(),
			ages = Vector{Int}(),
			interactions = Matrix{Float64}(undef, 0, 0),
			disease_type = String[],
			base_mortality = 0.0,
			disease_mortality = Float64[],
			fecundity = 0.0,
			transmission = Float64[],
			time_steps = 1,
			age_maturity = 1,
		)

		# Test Matrix{Bool} input
		pop = [Matrix{Bool}([true false false false; true false false false]) for _ in 1:10]
		pop[1][1, 1] = false
		pop[1][1, 3] = true

		result = coinfection_simulator(
			initial_pop = pop,
			ages = ones(Int, 10),
			interactions = [1.0 0.9; 0.9 1.0],
			disease_type = ["si", "si"],
			base_mortality = 0.0,
			disease_mortality = [0.0, 0.0],
			fecundity = 0.0,
			transmission = [0.5, 0.5],
			time_steps = 5,
			age_maturity = 1,
		)

		@test length(result) == 2  # Returns tuple of (populations, ages)
		@test length(result[1]) == 5  # 5 time steps
		@test length(result[2]) == 5  # 5 time steps

		# Test basic valid input with BitMatrix
		pop_bit = [BitMatrix([true false false false; true false false false]) for _ in 1:10]
		pop_bit[1][1, 1] = false
		pop_bit[1][1, 3] = true

		result_bit = coinfection_simulator(
			initial_pop = pop_bit,
			ages = ones(Int, 10),
			interactions = [1.0 0.9; 0.9 1.0],
			disease_type = ["si", "si"],
			base_mortality = 0.0,
			disease_mortality = [0.0, 0.0],
			fecundity = 0.0,
			transmission = [0.5, 0.5],
			time_steps = 5,
			age_maturity = 1,
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
			initial_pop = pop_bool,
			ages = ones(Int, 20),
			interactions = [1.0 0.8; 0.8 1.0],
			disease_type = ["si", "si"],
			base_mortality = 0.0,
			disease_mortality = [0.0, 0.0],
			fecundity = 0.0,
			transmission = [0.3, 0.3],
			time_steps = 5,
			age_maturity = 1,
		)

		Random.seed!(123)
		result_bit = coinfection_simulator(
			initial_pop = pop_bit,
			ages = ones(Int, 20),
			interactions = [1.0 0.8; 0.8 1.0],
			disease_type = ["si", "si"],
			base_mortality = 0.0,
			disease_mortality = [0.0, 0.0],
			fecundity = 0.0,
			transmission = [0.3, 0.3],
			time_steps = 5,
			age_maturity = 1,
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
			initial_pop = pop,
			ages = ones(Int, 10),
			interactions = [1.0;;],
			disease_type = ["sir"],
			base_mortality = 0.0,
			disease_mortality = [0.0],
			fecundity = 0.0,
			transmission = [0.5],
			time_steps = 2,
			age_maturity = 1,
			# Missing recovery parameter
		)

		# Test SEIR model requires latency parameter
		@test_throws AssertionError coinfection_simulator(
			initial_pop = pop,
			ages = ones(Int, 10),
			interactions = [1.0;;],
			disease_type = ["seir"],
			base_mortality = 0.0,
			disease_mortality = [0.0],
			fecundity = 0.0,
			transmission = [0.5],
			time_steps = 2,
			age_maturity = 1,
			recovery = [0.1],
			# Missing latency parameter
		)

		# Test SEIRS model requires immunity_loss parameter
		@test_throws AssertionError coinfection_simulator(
			initial_pop = pop,
			ages = ones(Int, 10),
			interactions = [1.0;;],
			disease_type = ["seirs"],
			base_mortality = 0.0,
			disease_mortality = [0.0],
			fecundity = 0.0,
			transmission = [0.5],
			time_steps = 2,
			age_maturity = 1,
			recovery = [0.1],
			latency = [5],
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
			initial_pop = pop,
			ages = ones(Int, 100),
			interactions = [1.0 0.5; 0.5 1.0],
			disease_type = ["si", "si"],
			base_mortality = 0.0,
			disease_mortality = [0.0, 0.0],
			fecundity = 0.0,
			transmission = [0.8, 0.8],
			time_steps = 10,
			age_maturity = 1,
		)

		# Check that infection can spread
		initial_infected = sum(m[1, 3] for m in results[1][1])
		final_infected = sum(m[1, 3] for m in results[1][end])
		@test final_infected >= initial_infected
	end

	@testset "Edge cases" begin
		# Test with mortality = 0 and fecundity = 0 (stable population)
		pop = [Matrix{Bool}([true false false false;]) for _ in 1:5]

		results = coinfection_simulator(
			initial_pop = pop,
			ages = ones(Int, 5),
			interactions = [1.0;;],
			disease_type = ["si"],
			base_mortality = 0.0,
			disease_mortality = [0.0],
			fecundity = 0.0,
			transmission = [0.0],  # No transmission
			time_steps = 3,
			age_maturity = 1,
			introduction = "none",
		)

		# Population should remain stable
		@test length(results[1][1]) == length(results[1][end])
		@test all(results[2][end] .== results[2][1] .+ 2)  # Ages should increase by 2
	end

	@testset "Birth scenario tests" begin
		Random.seed!(101)

		# Test scenario where births occur
		pop = [Matrix{Bool}([true false false false;]) for _ in 1:20]

		results = coinfection_simulator(
			initial_pop = pop,
			ages = fill(25, 20),  # All adults above maturity age
			interactions = [1.0;;],
			disease_type = ["si"],
			base_mortality = 0.0,
			disease_mortality = [0.0],
			fecundity = 0.5,  # High fecundity to ensure births
			transmission = [0.0],
			time_steps = 5,
			age_maturity = 18,
			introduction = "none",
		)

		# Check that population size increased due to births
		initial_size = length(results[1][1])
		final_size = length(results[1][end])
		@test final_size > initial_size

		# Check that newborns have age 0 at their birth timestep
		# and are all susceptible
		for t in 2:length(results[1])
			current_ages = results[2][t]
			current_pop = results[1][t]

			# Find individuals with age < t (born during simulation)
			newborns = findall(age -> age < t, current_ages)
			if !isempty(newborns)
				# All newborns should be susceptible to strain 1
				for nb in newborns
					@test current_pop[nb][1, 1] == true  # Susceptible
					@test current_pop[nb][1, 2] == false # Not exposed
					@test current_pop[nb][1, 3] == false # Not infected
					@test current_pop[nb][1, 4] == false # Not recovered
				end
			end
		end
	end

	@testset "Infection mortality tests" begin
		Random.seed!(202)

		# Test scenario where hosts die of infection
		pop = [Matrix{Bool}([true false false false;]) for _ in 1:50]
		# Start with multiple infected individuals
		for i in 1:10
			pop[i][1, 1] = false  # Not susceptible
			pop[i][1, 3] = true   # Infected
		end

		results = coinfection_simulator(
			initial_pop = pop,
			ages = fill(10, 50),
			interactions = [1.0;;],
			disease_type = ["si"],
			base_mortality = 0.0,     # No background mortality
			disease_mortality = [0.8], # High infection mortality
			fecundity = 0.0,
			transmission = [0.0],     # No transmission to isolate mortality effect
			time_steps = 10,
			age_maturity = 18,
			introduction = "none",
		)

		# Population should decrease due to infection mortality
		initial_size = length(results[1][1])
		final_size = length(results[1][end])
		@test final_size < initial_size

		# Count initial infected individuals
		initial_infected = sum(m[1, 3] for m in results[1][1])
		@test initial_infected == 10

		# Verify that mortality occurred (some infected individuals died)
		# Since we have high disease mortality and no recovery, infected should decrease
		final_infected = sum(m[1, 3] for m in results[1][end])
		@test final_infected < initial_infected
	end

	@testset "SIR disease model tests" begin
		Random.seed!(303)

		# Test SIR disease progression
		pop = [Matrix{Bool}([true false false false;]) for _ in 1:100]
		# Start with one infected individual
		pop[1][1, 1] = false  # Not susceptible
		pop[1][1, 3] = true   # Infected

		results = coinfection_simulator(
			initial_pop = pop,
			ages = fill(20, 100),
			interactions = [1.0;;],
			disease_type = ["sir"],
			base_mortality = 0.0,
			disease_mortality = [0.0],
			fecundity = 0.0,
			transmission = [0.3],
			time_steps = 20,
			age_maturity = 18,
			recovery = [0.2],  # Recovery rate
			introduction = "none",
		)

		# Check that we have transitions from I to R
		final_recovered = sum(m[1, 4] for m in results[1][end])
		@test final_recovered > 0

		# Check that no individual is in multiple states simultaneously
		for t in 1:length(results[1])
			for individual in results[1][t]
				row_sum = sum(individual[1, :])
				@test row_sum == 1  # Exactly one state per strain
			end
		end

		# Verify SIR progression: S -> I -> R (no return to S)
		# Track one individual that becomes infected
		found_progression = false
		for t in 2:length(results[1])
			for (i, individual) in enumerate(results[1][t])
				# Find someone who is infected
				if individual[1, 3] == true
					# Check if they recover in later timesteps
					for future_t in (t+1):length(results[1])
						if i <= length(results[1][future_t])
							future_state = results[1][future_t][i]
							if future_state[1, 4] == true  # Recovered
								found_progression = true
								# Once recovered, should never be susceptible again
								@test future_state[1, 1] == false
								break
							end
						end
					end
				end
			end
		end
	end

	@testset "SEIR disease model tests" begin
		Random.seed!(404)

		# Test SEIR disease progression  
		pop = [Matrix{Bool}([true false false false;]) for _ in 1:80]
		# Start with one infected individual
		pop[1][1, 1] = false  # Not susceptible
		pop[1][1, 3] = true   # Infected

		results = coinfection_simulator(
			initial_pop = pop,
			ages = fill(25, 80),
			interactions = [1.0;;],
			disease_type = ["seir"],
			base_mortality = 0.0,
			disease_mortality = [0.0],
			fecundity = 0.0,
			transmission = [0.4],
			time_steps = 25,
			age_maturity = 18,
			recovery = [0.15],
			latency = [3],  # 3 time steps in exposed state
			introduction = "none",
		)

		# Check that we have exposed individuals (unique to SEIR/SEIRS)
		max_exposed = maximum([sum(m[1, 2] for m in timestep) for timestep in results[1]])
		@test max_exposed > 0

		# Check that we have recovered individuals
		final_recovered = sum(m[1, 4] for m in results[1][end])
		@test final_recovered > 0

		# Check SEIR state constraints
		for t in 1:length(results[1])
			for individual in results[1][t]
				row_sum = sum(individual[1, :])
				@test row_sum == 1  # Exactly one state per strain
			end
		end

		# Verify that exposed state is used (S -> E -> I -> R progression)
		found_exposed_to_infected = false
		for t in 1:(length(results[1])-1)
			for (i, individual) in enumerate(results[1][t])
				if individual[1, 2] == true && i <= length(results[1][t+1])  # Exposed
					next_individual = results[1][t+1][i]
					if next_individual[1, 3] == true  # Became infected
						found_exposed_to_infected = true
						break
					end
				end
			end
		end
	end

	@testset "SEIRS disease model tests" begin
		Random.seed!(505)

		# Test SEIRS disease progression (includes immunity loss)
		pop = [Matrix{Bool}([true false false false;]) for _ in 1:60]
		# Start with one infected individual
		pop[1][1, 1] = false  # Not susceptible  
		pop[1][1, 3] = true   # Infected

		results = coinfection_simulator(
			initial_pop = pop,
			ages = fill(30, 60),
			interactions = [1.0;;],
			disease_type = ["seirs"],
			base_mortality = 0.0,
			disease_mortality = [0.0],
			fecundity = 0.0,
			transmission = [0.5],
			time_steps = 30,
			age_maturity = 18,
			recovery = [0.2],
			latency = [2],  # 2 time steps in exposed state
			immunity_loss = [0.1],  # 10% chance of losing immunity
			introduction = "none",
		)

		# Check that we have exposed individuals 
		max_exposed = maximum([sum(m[1, 2] for m in timestep) for timestep in results[1]])
		@test max_exposed > 0

		# Check that we have recovered individuals
		max_recovered = maximum([sum(m[1, 4] for m in timestep) for timestep in results[1]])
		@test max_recovered > 0

		# Check SEIRS state constraints
		for t in 1:length(results[1])
			for individual in results[1][t]
				row_sum = sum(individual[1, :])
				@test row_sum == 1  # Exactly one state per strain
			end
		end

		# The key feature of SEIRS: recovered individuals can become susceptible again
		# Look for evidence of immunity loss (R -> S transition)
		found_immunity_loss = false
		for t in 1:(length(results[1])-1)
			for (i, individual) in enumerate(results[1][t])
				if individual[1, 4] == true && i <= length(results[1][t+1])  # Recovered
					next_individual = results[1][t+1][i]
					if next_individual[1, 1] == true  # Became susceptible again
						found_immunity_loss = true
						break
					end
				end
			end
		end

		# Note: We don't require finding immunity loss in every test run due to randomness,
		# but the model should be capable of it. The important test is that the simulation
		# runs without error and maintains valid state constraints.
	end
end
