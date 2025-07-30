# Advanced Example: Complex Disease Models and Strain Interactions
#
# This example demonstrates advanced features of CoinfectionSimulator.jl:
# - Different disease models (SI, SIR, SEIR, SEIRS)
# - Complex interaction matrices with facilitation and competition
# - Strain introduction scheduling
# - Demographic processes (births, deaths, aging)
# - Detailed analysis and comparison

using CoinfectionSimulator
using Random
using Plots
using DataFrames
using Statistics

println("=== Advanced Multi-Strain Coinfection Simulation ===\n")

# Set seed for reproducibility
Random.seed!(123)

## 1. Define Multiple Disease Models
println("Setting up diverse disease models...")

# Strain 1: Fast-spreading SI model (like a persistent infection)
strain1_model = SIModel(
    0.25,    # High transmission rate
    0.001    # Low mortality
)

# Strain 2: Classic SIR model (like influenza)
strain2_model = SIRModel(
    0.20,    # Moderate transmission
    0.003,   # Low mortality
    0.15     # Recovery rate
)

# Strain 3: SEIR model with latency (like COVID-19)
strain3_model = SEIRModel(
    0.18,    # Transmission rate
    0.005,   # Mortality rate
    0.12,    # Recovery rate
    5        # 5-day latency period
)

# Strain 4: SEIRS model with immunity loss (like cold virus)
strain4_model = SEIRSModel(
    0.15,    # Transmission rate
    0.001,   # Very low mortality
    0.20,    # Fast recovery
    3,       # Short latency
    0.05     # Immunity loss rate
)

disease_models = [strain1_model, strain2_model, strain3_model, strain4_model]
println("✓ Created 4 different disease models:")
println("  - Strain 1: SI (persistent)")
println("  - Strain 2: SIR (acute with immunity)")
println("  - Strain 3: SEIR (with incubation)")
println("  - Strain 4: SEIRS (immunity wanes)")

## 2. Create Complex Interaction Matrix
println("\nCreating complex strain interaction matrix...")

# Create a matrix demonstrating various interaction types
interactions = create_interaction_matrix(
    4,        # 4 strains
    true,     # Asymmetric (priority effects)
    0.4,      # Moderate interaction strength
    cf_ratio=0.6  # More facilitation than competition
)

println("✓ Generated interaction matrix:")
display(round.(interactions, digits=3))

# Interpret the matrix
println("\nInteraction interpretation:")
for i in 1:4, j in 1:4
    if i != j
        effect = interactions[i, j]
        if effect > 1.0
            println("  Strain $i FACILITATES strain $j (×$(round(effect, digits=2)))")
        elseif effect < 1.0
            println("  Strain $i COMPETES with strain $j (×$(round(effect, digits=2)))")
        end
    end
end

## 3. Initialize Large Population
println("\nInitializing population...")

n_individuals = 500
population = Population(Individual[])

# Create age-structured population
for i in 1:n_individuals
    # Realistic age distribution (more young people)
    age = if rand() < 0.3
        rand(0:18)     # 30% children/teens
    elseif rand() < 0.5
        rand(19:35)    # 35% young adults
    elseif rand() < 0.3
        rand(36:60)    # 30% middle-aged
    else
        rand(61:85)    # 15% elderly
    end

    individual = Individual(4, age)  # 4 strains, all susceptible
    push!(population.individuals, individual)
end

# Initial seeding - only strain 1 present initially
for i in 1:3
    population[i][1, 1] = false  # No longer susceptible
    population[i][1, 3] = true   # Infected with strain 1
end

println("✓ Created population of $n_individuals individuals")
println("  - Age distribution: realistic (more young people)")
println("  - Initial infections: 3 individuals with strain 1")

## 4. Set Advanced Simulation Parameters
println("\nConfiguring advanced simulation parameters...")

params = SimulationParameters(
    disease_models,
    interactions,
    0.0005,              # Very low base mortality
    0.025,               # Birth rate
    18,                  # Age of maturity
    :random,             # Random strain introductions
    200                  # Longer simulation
)

println("✓ Advanced parameters set:")
println("  - Duration: 200 time steps")
println("  - Demographic processes: births and deaths")
println("  - Random strain introductions enabled")

## 5. Run Extended Simulation
println("\nRunning extended simulation...")

results = simulate(population, params)
populations, age_vectors = results

println("✓ Simulation completed!")
println("  - Initial population: $(length(populations[1]))")
println("  - Final population: $(length(populations[end]))")

## 6. Comprehensive Analysis
println("\nPerforming comprehensive analysis...")

n_timesteps = length(populations)
n_strains = 4

# Initialize tracking arrays
susceptible = zeros(Int, n_timesteps, n_strains)
exposed = zeros(Int, n_timesteps, n_strains)
infected = zeros(Int, n_timesteps, n_strains)
recovered = zeros(Int, n_timesteps, n_strains)
total_pop = zeros(Int, n_timesteps)
coinfections = zeros(Int, n_timesteps, n_strains, n_strains)

# Track all states over time
for t in 1:n_timesteps
    pop = populations[t]
    total_pop[t] = length(pop)

    for individual in pop
        for strain in 1:n_strains
            if individual[strain, 1]
                susceptible[t, strain] += 1
            elseif individual[strain, 2]
                exposed[t, strain] += 1
            elseif individual[strain, 3]
                infected[t, strain] += 1
            elseif individual[strain, 4]
                recovered[t, strain] += 1
            end
        end

        # Track coinfections
        for i in 1:n_strains, j in (i+1):n_strains
            if individual[i, 3] && individual[j, 3]
                coinfections[t, i, j] += 1
            end
        end
    end
end

# Calculate epidemic metrics
println("\n=== EPIDEMIC SUMMARY ===")
for strain in 1:n_strains
    peak_infected = maximum(infected[:, strain])
    peak_time = findfirst(infected[:, strain] .== peak_infected)
    final_recovered = recovered[end, strain]
    attack_rate = final_recovered / n_individuals * 100

    println("Strain $strain:")
    println("  - Peak infections: $peak_infected (at time $peak_time)")
    println("  - Final recovered: $final_recovered")
    println("  - Attack rate: $(round(attack_rate, digits=1))%")
end

# Coinfection analysis
println("\nCoinfection patterns:")
for i in 1:n_strains, j in (i+1):n_strains
    max_coinfection = maximum(coinfections[:, i, j])
    if max_coinfection > 0
        peak_time = findfirst(coinfections[:, i, j] .== max_coinfection)
        println("  - Strains $i & $j: peak $max_coinfection (at time $peak_time)")
    end
end

## 7. Create Advanced Visualizations
println("\nCreating advanced visualizations...")

time_steps = 1:n_timesteps

# Plot 1: All strain infections over time
p1 = plot(title="Multi-Strain Infection Dynamics", xlabel="Time Step",
    ylabel="Infected Individuals", legend=:topright)
colors = [:red, :blue, :green, :orange]
for strain in 1:n_strains
    plot!(p1, time_steps, infected[:, strain],
        label="Strain $strain", lw=2, color=colors[strain])
end

# Plot 2: Disease states for most successful strain
dominant_strain = argmax([maximum(infected[:, i]) for i in 1:n_strains])
p2 = plot(title="Disease States - Strain $dominant_strain (Dominant)",
    xlabel="Time Step", ylabel="Number of Individuals")
plot!(p2, time_steps, susceptible[:, dominant_strain],
    label="Susceptible", lw=2, color=:lightgray)
plot!(p2, time_steps, exposed[:, dominant_strain],
    label="Exposed", lw=2, color=:yellow)
plot!(p2, time_steps, infected[:, dominant_strain],
    label="Infected", lw=2, color=:red)
plot!(p2, time_steps, recovered[:, dominant_strain],
    label="Recovered", lw=2, color=:darkgreen)

# Plot 3: Population demographics
p3 = plot(title="Population Dynamics", xlabel="Time Step",
    ylabel="Number of Individuals")
plot!(p3, time_steps, total_pop, label="Total Population",
    lw=3, color=:black)

# Add birth and death rates (approximate)
births = diff(total_pop)
births[births.<0] .= 0  # Only count net increases as births
deaths = -diff(total_pop)
deaths[deaths.<0] .= 0  # Only count net decreases as deaths

if length(births) > 0
    plot!(p3, time_steps[2:end], cumsum(births),
        label="Cumulative Births", lw=2, color=:green, linestyle=:dash)
    plot!(p3, time_steps[2:end], cumsum(deaths),
        label="Cumulative Deaths", lw=2, color=:purple, linestyle=:dash)
end

# Plot 4: Coinfection heatmap (average over time)
avg_coinfections = zeros(n_strains, n_strains)
for i in 1:n_strains, j in 1:n_strains
    if i != j
        avg_coinfections[i, j] = mean(coinfections[:, min(i, j), max(i, j)])
    end
end

p4 = heatmap(avg_coinfections,
    title="Average Coinfection Levels",
    xlabel="Strain", ylabel="Strain",
    color=:viridis)

# Combine all plots
final_plot = plot(p1, p2, p3, p4,
    layout=(2, 2),
    size=(1200, 1000),
    plot_title="Advanced Multi-Strain Simulation Analysis")

savefig(final_plot, "advanced_coinfection_analysis.png")
println("✓ Saved comprehensive analysis as 'advanced_coinfection_analysis.png'")

## 8. Comparative Analysis with Different Interaction Strengths
println("\nRunning comparative analysis...")

# Compare with no interactions
interactions_neutral = ones(4, 4)
params_neutral = SimulationParameters(
    disease_models, interactions_neutral, 0.0005, 0.025, 18, :random, 200
)

Random.seed!(123)  # Same seed for comparison
results_neutral = simulate(deepcopy(population), params_neutral)

# Quick comparison
infected_neutral = zeros(Int, length(results_neutral[1]), n_strains)
for t in 1:length(results_neutral[1])
    pop = results_neutral[1][t]
    for individual in pop, strain in 1:n_strains
        if individual[strain, 3]
            infected_neutral[t, strain] += 1
        end
    end
end

println("\nComparison: Complex Interactions vs No Interactions")
println("Peak infections with complex interactions:")
for strain in 1:n_strains
    peak_complex = maximum(infected[:, strain])
    peak_neutral = maximum(infected_neutral[:, strain])
    change = peak_complex - peak_neutral
    pct_change = round((change / peak_neutral) * 100, digits=1)
    println("  Strain $strain: $peak_complex vs $peak_neutral ($(change >= 0 ? "+" : "")$change, $(change >= 0 ? "+" : "")$pct_change%)")
end

## 9. Age-Stratified Analysis
println("\nAge-stratified analysis...")

# Analyze final outcomes by age group
age_groups = ["0-18", "19-35", "36-60", "61+"]
age_ranges = [(0, 18), (19, 35), (36, 60), (61, 100)]

final_pop = populations[end]
age_analysis = DataFrame(
    AgeGroup=String[],
    Count=Int[],
    AnyInfection=Int[],
    Strain1Ever=Int[],
    Strain2Ever=Int[],
    Strain3Ever=Int[],
    Strain4Ever=Int[]
)

for (i, (min_age, max_age)) in enumerate(age_ranges)
    age_individuals = filter(ind -> min_age <= ind.age <= max_age, final_pop)
    count = length(age_individuals)

    any_infection = sum(any(ind[strain, 3] || ind[strain, 4] for strain in 1:4) for ind in age_individuals)
    strain_infections = [sum(ind[strain, 3] || ind[strain, 4] for ind in age_individuals) for strain in 1:4]

    push!(age_analysis, (
        age_groups[i], count, any_infection,
        strain_infections[1], strain_infections[2],
        strain_infections[3], strain_infections[4]
    ))
end

println("\nAge-stratified infection rates:")
display(age_analysis)

## 10. Export Results
println("\nExporting detailed results...")

# Create detailed time series data
results_df = DataFrame(
    TimeStep=Int[],
    TotalPopulation=Int[],
    Strain1_S=Int[], Strain1_E=Int[], Strain1_I=Int[], Strain1_R=Int[],
    Strain2_S=Int[], Strain2_E=Int[], Strain2_I=Int[], Strain2_R=Int[],
    Strain3_S=Int[], Strain3_E=Int[], Strain3_I=Int[], Strain3_R=Int[],
    Strain4_S=Int[], Strain4_E=Int[], Strain4_I=Int[], Strain4_R=Int[]
)

for t in 1:n_timesteps
    push!(results_df, (
        t, total_pop[t],
        susceptible[t, 1], exposed[t, 1], infected[t, 1], recovered[t, 1],
        susceptible[t, 2], exposed[t, 2], infected[t, 2], recovered[t, 2],
        susceptible[t, 3], exposed[t, 3], infected[t, 3], recovered[t, 3],
        susceptible[t, 4], exposed[t, 4], infected[t, 4], recovered[t, 4]
    ))
end

# Save to CSV
using CSV
CSV.write("advanced_simulation_results.csv", results_df)
println("✓ Detailed results saved to 'advanced_simulation_results.csv'")

println("\n=== Advanced Example Complete! ===")
println("This example demonstrated:")
println("• Multiple disease models (SI, SIR, SEIR, SEIRS)")
println("• Complex asymmetric interaction matrices")
println("• Age-structured populations")
println("• Demographic processes (births, deaths)")
println("• Random strain introductions")
println("• Comprehensive epidemic analysis")
println("• Coinfection pattern tracking")
println("• Comparative analysis")
println("• Age-stratified outcomes")
println("• Data export capabilities")
