# Basic Example: CoinfectionSimulator.jl
#
# This example demonstrates the core functionality of CoinfectionSimulator.jl,
# showing how to set up and run a basic two-strain coinfection simulation
# using the modern API.

using CoinfectionSimulator
using Random
using Plots

# Set random seed for reproducibility
Random.seed!(42)

## 1. Define Disease Models

# Create SIR models for two different strains
strain1_model = SIRModel(
    0.15,    # transmission rate
    0.005,   # disease mortality
    0.1      # recovery rate
)

strain2_model = SIRModel(
    0.12,    # slightly lower transmission
    0.003,   # lower disease mortality
    0.08     # slower recovery
)

disease_models = [strain1_model, strain2_model]

## 2. Create Interaction Matrix

# Create a matrix with slight competition between strains
# interactions[i,j] = effect of strain i on transmission of strain j
interactions = [1.0 0.9;   # strain 2 slightly reduces strain 1 transmission
    0.9 1.0]   # strain 1 slightly reduces strain 2 transmission

println("✓ Interaction matrix:")
display(interactions)

## 3. Initialize Population
println("\nInitializing population...")

n_individuals = 200
population = Population(Individual[])

# Create individuals - all start susceptible to both strains
for i in 1:n_individuals
    age = rand(1:80)  # Random age between 1-80
    individual = Individual(2, age)  # 2 strains, all susceptible initially
    push!(population.individuals, individual)
end

# Seed the epidemic by infecting a few individuals with strain 1
for i in 1:5
    population[i][1, 1] = false  # No longer susceptible to strain 1
    population[i][1, 3] = true   # Now infected with strain 1
end

# Seed strain 2 as well (introduced later)
population[n_individuals][2, 1] = false  # No longer susceptible to strain 2
population[n_individuals][2, 3] = true   # Now infected with strain 2

## 4. Set Simulation Parameters

params = SimulationParameters(
    disease_models,           # Our SIR models
    interactions,            # Strain interactions
    0.001,                   # Base mortality rate
    0.02,                    # Fecundity (birth rate)
    18,                      # Age of maturity
    :none,                   # No additional introductions
    150                      # Number of time steps
)

## 5. Run Simulation

results = simulate(population, params)
populations, age_vectors = results

## 6. Analyze Results

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

    for individual in pop
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

println("=== Simulation Summary ===")
println("Peak infections:")
println("  - Strain 1: $peak_infections_1")
println("  - Strain 2: $peak_infections_2")
println("  - Coinfections: $peak_coinfections")
println("Final recovered:")
println("  - Strain 1: $final_recovered_1")
println("  - Strain 2: $final_recovered_2")
println("Population change: $(length(populations[1])) → $(length(populations[end]))")

## 7. Create Visualizations

time_steps = 1:n_timesteps

# Plot 1: Infection dynamics
p1 = plot(time_steps, [infected_1, infected_2, coinfected],
    label=["Strain 1 Infected" "Strain 2 Infected" "Coinfected"],
    title="Infection Dynamics",
    xlabel="Time Step",
    ylabel="Number of Individuals",
    lw=2,
    color=[:red :blue :purple])

# Plot 2: Recovery dynamics
p2 = plot(time_steps, [recovered_1, recovered_2],
    label=["Strain 1 Recovered" "Strain 2 Recovered"],
    title="Recovery Dynamics",
    xlabel="Time Step",
    ylabel="Number of Individuals",
    lw=2,
    color=[:darkred :darkblue])

# Plot 3: Population size and susceptible individuals
p3 = plot(time_steps, total_pop,
    label="Total Population",
    title="Population Dynamics",
    xlabel="Time Step",
    ylabel="Number of Individuals",
    lw=2,
    color=:black)

plot!(p3, time_steps, susceptible_1,
    label="Susceptible to Strain 1",
    lw=2,
    color=:lightcoral,
    linestyle=:dash)

plot!(p3, time_steps, susceptible_2,
    label="Susceptible to Strain 2",
    lw=2,
    color=:lightblue,
    linestyle=:dash)

# Combine plots
final_plot = plot(p1, p2, p3,
    layout=(3, 1),
    size=(800, 900),
    plot_title="Two-Strain Coinfection Simulation")

# Save visualization
savefig(final_plot, "basic_coinfection_example.png")

## 8. Demonstrate Virtual Ecologist Sampling

# Sample the population over time to simulate surveillance
sampling_params = SamplingParameters(
    0.1,    # Sample 10% of population
    0.05,   # 5% false positive rate
    0.1     # 10% false negative rate
)

detection_results = sample_populations(populations, sampling_params)

# Count time steps where each strain was detected
strain1_detections = sum(detection_results[:, 1])
strain2_detections = sum(detection_results[:, 2])

println("Sampling results (10% of population sampled):")
println("  - Strain 1 detected in $strain1_detections/$(n_timesteps) time steps")
println("  - Strain 2 detected in $strain2_detections/$(n_timesteps) time steps")

# Show first 20 time steps of detection
println("\nFirst 20 time steps - Detection pattern:")
println("Time\tStrain1\tStrain2")
for t in 1:min(20, n_timesteps)
    s1 = detection_results[t, 1] ? "✓" : "✗"
    s2 = detection_results[t, 2] ? "✓" : "✗"
    println("$t\t$s1\t$s2")
end
