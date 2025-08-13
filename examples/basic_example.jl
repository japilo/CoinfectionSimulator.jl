# # [Basic Example: CoinfectionSimulator.jl](@id basic_example)
#
# This example demonstrates the core functionality of CoinfectionSimulator.jl,
# showing how to set up and run a basic two-strain coinfection simulation with
# slight competition between the strains.
#
# ## Biological Context
# This simulation models two co-circulating strains that compete for hosts.
# The competition could result from:
# - Cross-reactive immune responses between related pathogens
# - Resource competition within infected hosts
# - Behavioral changes to the host that affect transmission patterns

using CoinfectionSimulator
using Random
using Plots
Random.seed!(42)

# ## 1. Define Disease Models

# CoinfectionSimulator.jl supports multiple compartmental disease models:
# - **SI**: Susceptible → Infected (no recovery, chronic infections)
# - **SIR**: Susceptible → Infected → Recovered (acute with lifelong immunity)
# - **SEIR**: Susceptible → Exposed → Infected → Recovered (includes latent period)
# - **SEIRS**: Susceptible → Exposed → Infected → Recovered → Susceptible (temporary immunity)
#
# For this example, we use SIR models for both strains, representing acute infections
# like influenza where individuals recover and develop strain-specific immunity.

# **Strain 1: More aggressive pathogen**
# Higher transmission but also higher mortality
strain1_model = SIRModel(
	0.15,    # transmission rate: 15% chance per infected contact per time step
	0.005,   # disease-induced mortality: 0.5% death rate for infected individuals
	0.1,     # recovery rate: average infectious period = 1/0.1 = 10 time steps
)

# **Strain 2: Milder but persistent pathogen**
# Lower transmission and mortality but longer infectious period
strain2_model = SIRModel(
	0.12,    # transmission rate: 12% per contact (slightly less transmissible)
	0.003,   # disease-induced mortality: 0.3% death rate (less virulent)
	0.08,    # recovery rate: average infectious period = 1/0.08 = 12.5 time steps
)

disease_models = [strain1_model, strain2_model]

# ## 2. Create Interaction Matrix

# The interaction matrix defines how the presence of one strain affects the transmission
# of another strain in potential new hosts. Each element [i,j] represents how strain j
# affects the transmission rate of strain i:
# - 1.0 = No interaction (baseline transmission)
# - < 1.0 = Competition (reduced transmission)
# - > 1.0 = Facilitation (enhanced transmission)
#
# **Biological Mechanisms:**
# - Competition: Cross-reactive immunity, resource competition within the host
# - Facilitation: Immunosuppression of the host, behavioral changes that facilitate further infections
#
# This example models **symmetric competition** where each strain reduces the other's
# transmission by 10%. The symmetry indicates that there are no priority effects.
interactions = [1.0 0.9;   # Row 1: strain 2 reduces strain 1 transmission by 10%
	            0.9 1.0]   # Row 2: strain 1 reduces strain 2 transmission by 10%

# ## 3. Initialize Population
# The package has a `Population` type that holds all individuals in the simulation. We need to define a population,
# including the number of individuals, their ages (in whatever unit the timestep of the simulation is meant to be)
# and their initial states. Each individual is represented by an `Individual` type, which are susceptible to all strains
# by default.

n_individuals = 200
population = Population(Individual[])

# **Create age-structured population**
# All individuals start susceptible to both strains
for i in 1:n_individuals
	age = rand(1:80)  # Random age 1-80 time steps (could represent days/weeks/months)
	individual = Individual(2, age)  # 2 strains, starts susceptible to both
	push!(population.individuals, individual)
end

# **Seed initial infections to start the epidemic**

# Strain 1 introduction (index cases)
for i in 1:5
	population[i][1, 1] = false  # No longer susceptible to strain 1
	population[i][1, 3] = true   # Now infected with strain 1
end

# Strain 2 introduction
population[n_individuals][2, 1] = false  # No longer susceptible to strain 2
population[n_individuals][2, 3] = true   # Now infected with strain 2

# ## 4. Set Simulation Parameters
# The `SimulationParameters` type holds all the parameters needed to run the simulation.

params = SimulationParameters(
	disease_models,           # Our SIR models (one per strain)
	interactions,            # Strain interactions
	0.001,                   # Base mortality rate
	0.02,                    # Fecundity (birth rate, supplied as lambda to a Poisson process)
	30,                      # Age of maturity
	:none,                   # No additional introductions of strains
	150,                      # Number of time steps
)

## 5. Run Simulation

# The simulate() function returns a tuple containing:
# - populations: Vector of Population objects (one per time step)
# - age_vectors: Vector of age distributions over time
#
# **What happens during simulation:**
# 1. Disease transmission from infected to susceptible individuals
# 2. Disease progression (recovery, mortality) according to model parameters
# 3. Demographic processes (births, deaths, aging)
# 4. Strain interactions modify transmission rates in coinfected hosts

populations = simulate(population, params)

# ## 6. Analyze Results
# The simulation returns a list of populations at each time step and their age vectors.
# Each population is a list of individuals, where each individual has a state vector for each strain
# (e.g., susceptible, exposed, infected, recovered).
# We will gather this information to analyze and visualize the dynamics of the two strains.

n_timesteps = length(populations)

# Track disease states over time
susceptible_1 = zeros(Int, n_timesteps)
exposed_1 = zeros(Int, n_timesteps)
infected_1 = zeros(Int, n_timesteps)
recovered_1 = zeros(Int, n_timesteps)

susceptible_2 = zeros(Int, n_timesteps)
exposed_2 = zeros(Int, n_timesteps)
infected_2 = zeros(Int, n_timesteps)
recovered_2 = zeros(Int, n_timesteps)

coinfected = zeros(Int, n_timesteps)
total_pop = zeros(Int, n_timesteps)

for t in 1:n_timesteps
	pop = populations[t]
	total_pop[t] = length(pop)

	for i in 1:length(pop.individuals)
		individual = pop.individuals[i]
		# Strain 1 states
		if individual[1, 1]
			susceptible_1[t] += 1
		end
		if individual[1, 2]
			exposed_1[t] += 1
		end
		if individual[1, 3]
			infected_1[t] += 1
		end
		if individual[1, 4]
			recovered_1[t] += 1
		end

		# Strain 2 states
		if individual[2, 1]
			susceptible_2[t] += 1
		end
		if individual[2, 2]
			exposed_2[t] += 1
		end
		if individual[2, 3]
			infected_2[t] += 1
		end
		if individual[2, 4]
			recovered_2[t] += 1
		end

		# Coinfection (infected with both strains)
		if individual[1, 3] && individual[2, 3]
			coinfected[t] += 1
		end
	end
end

# Calculate summary statistics
peak_infections_1 = maximum(infected_1)
peak_infections_2 = maximum(infected_2)
peak_coinfections = maximum(coinfected)
final_recovered_1 = recovered_1[end]
final_recovered_2 = recovered_2[end]

time_steps = 1:n_timesteps

# Plot 1: Infection dynamics
p1 = plot(time_steps, [infected_1, infected_2, coinfected],
	label = ["Strain 1 Infected" "Strain 2 Infected" "Coinfected"],
	title = "Infection Dynamics",
	xlabel = "Time Step",
	ylabel = "Number of Individuals",
	lw = 2,
	color = [:red :blue :purple])

# Plot 2: Recovery dynamics
p2 = plot(time_steps, [recovered_1, recovered_2],
	label = ["Strain 1 Recovered" "Strain 2 Recovered"],
	title = "Recovery Dynamics",
	xlabel = "Time Step",
	ylabel = "Number of Individuals",
	lw = 2,
	color = [:darkred :darkblue])

# Plot 3: Population size and susceptible individuals
p3 = plot(time_steps, total_pop,
	label = "Total Population",
	title = "Population Dynamics",
	xlabel = "Time Step",
	ylabel = "Number of Individuals",
	lw = 2,
	color = :black)

plot!(p3, time_steps, susceptible_1,
	label = "Susceptible to Strain 1",
	lw = 2,
	color = :lightcoral,
	linestyle = :dash)

plot!(p3, time_steps, susceptible_2,
	label = "Susceptible to Strain 2",
	lw = 2,
	color = :lightblue,
	linestyle = :dash)

# Combine plots
final_plot = plot(p1, p2, p3,
	layout = (3, 1),
	size = (800, 900),
	plot_title = "Two-Strain Coinfection Simulation")

# Save visualization
savefig(final_plot, "basic_coinfection_example.png")

# ## 8. Virtual Ecologist Sampling
# `CoinfectionSimulator.jl` also includes functionality for virtual ecologist sampling,
# which allows you to simulate imperfect sampling of the population over time to detect infections.

# Sample the population over time to simulate surveillance
sampling_params = SamplingParameters(
	0.1,    # Sample 10% of population
	0.05,   # 5% false positive rate
	0.1,     # 10% false negative rate
)

detection_results = sample_populations(populations, sampling_params)
