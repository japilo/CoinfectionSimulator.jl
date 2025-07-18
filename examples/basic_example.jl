# Basic Example: Two-Strain Coinfection Simulation

using CoinfectionSimulator
using Random
using Plots

# Set random seed for reproducibility
Random.seed!(42)

println("Setting up two-strain coinfection simulation...")

## Population Setup
n_individuals = 200
n_strains = 2

# Initialize population - everyone starts susceptible
initial_pop = [falses(n_strains, 4) for _ in 1:n_individuals]
for individual in initial_pop
    individual[:, 1] .= true  # All individuals start susceptible
end

# Random ages for the population
ages = rand(1:80, n_individuals)

## Simulation Parameters
# Interaction matrix: slight competition between strains
interactions = [1.0 0.9; 0.9 1.0]

# Disease parameters
disease_types = ["sir", "sir"]  # Both strains use SIR model
transmission = [0.15, 0.12]     # Different transmission rates
recovery = [0.1, 0.08]          # Different recovery rates

# Population parameters
base_mortality = 0.001          # Low background mortality
disease_mortality = [0.005, 0.003]  # Additional disease mortality
fecundity = 0.02               # Birth rate
age_maturity = 20              # Reproductive age

## Run Simulation
println("Running simulation for 100 time steps...")

results = coinfection_simulator(
    initial_pop=initial_pop,
    ages=ages,
    interactions=interactions,
    disease_type=disease_types,
    base_mortality=base_mortality,
    disease_mortality=disease_mortality,
    fecundity=fecundity,
    transmission=transmission,
    time_steps=100,
    age_maturity=age_maturity,
    recovery=recovery,
    introduction="simultaneous"
)

populations, age_vectors = results

## Analyze Results
println("Analyzing simulation results...")

# Count disease states over time
n_timesteps = length(populations)
susceptible_1 = zeros(Int, n_timesteps)
infected_1 = zeros(Int, n_timesteps)
recovered_1 = zeros(Int, n_timesteps)
susceptible_2 = zeros(Int, n_timesteps)
infected_2 = zeros(Int, n_timesteps)
recovered_2 = zeros(Int, n_timesteps)
total_pop = zeros(Int, n_timesteps)

for t in 1:n_timesteps
    pop = populations[t]
    total_pop[t] = length(pop)
    
    for individual in pop
        # Strain 1
        if individual[1, 1]; susceptible_1[t] += 1; end
        if individual[1, 3]; infected_1[t] += 1; end
        if individual[1, 4]; recovered_1[t] += 1; end
        
        # Strain 2  
        if individual[2, 1]; susceptible_2[t] += 1; end
        if individual[2, 3]; infected_2[t] += 1; end
        if individual[2, 4]; recovered_2[t] += 1; end
    end
end

# Plot results
time_steps = 1:n_timesteps

p1 = plot(time_steps, [infected_1, infected_2], 
         label=["Strain 1 Infected" "Strain 2 Infected"],
         title="Infection Dynamics", 
         xlabel="Time Step", 
         ylabel="Number of Individuals",
         lw=2)

p2 = plot(time_steps, total_pop,
         label="Total Population",
         title="Population Size",
         xlabel="Time Step",
         ylabel="Number of Individuals",
         lw=2, color=:black)

p3 = plot(time_steps, [recovered_1, recovered_2],
         label=["Strain 1 Recovered" "Strain 2 Recovered"],
         title="Recovery Dynamics",
         xlabel="Time Step",
         ylabel="Number of Individuals",
         lw=2)

plot(p1, p2, p3, layout=(3,1), size=(800, 600))
savefig("coinfection_dynamics.png")

println("Simulation completed!")
println("Final population size: $(total_pop[end])")
println("Peak infections - Strain 1: $(maximum(infected_1)), Strain 2: $(maximum(infected_2))")
println("Final recovered - Strain 1: $(recovered_1[end]), Strain 2: $(recovered_2[end])")

println("\nExample completed! Check 'coinfection_dynamics.png' for visualizations.")
