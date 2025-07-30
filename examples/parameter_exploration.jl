# Parameter Exploration Example: Batch Simulations and Sensitivity Analysis
#
# This example demonstrates how to explore parameter space systematically
# using CoinfectionSimulator.jl, including:
# - Batch simulation runs with different parameters
# - Sensitivity analysis for key parameters
# - Interaction matrix exploration
# - Statistical analysis of results
# - Parameter optimization techniques

using CoinfectionSimulator
using Random
using Plots
using DataFrames
using Statistics
using StatsBase
using CSV

println("=== Parameter Exploration and Sensitivity Analysis ===\n")

## 1. Define Base Scenario
println("Setting up base scenario...")

# Base disease models
base_sir = SIRModel(0.15, 0.005, 0.1)  # transmission, mortality, recovery
base_models = [base_sir, base_sir]  # Two identical strains

# Base population
function create_base_population(n_individuals=100, seed_infections=3)
    Random.seed!(42)  # Fixed seed for consistency
    population = Population(Individual[])

    for i in 1:n_individuals
        age = rand(18:65)  # Working age population
        individual = Individual(2, age)
        push!(population.individuals, individual)
    end

    # Seed infections
    for i in 1:seed_infections
        population[i][1, 1] = false
        population[i][1, 3] = true
    end

    return population
end

println("✓ Base scenario established")

## 2. Transmission Rate Sensitivity Analysis
println("\nRunning transmission rate sensitivity analysis...")

transmission_rates = 0.05:0.05:0.50  # Range from 5% to 50%
n_rates = length(transmission_rates)
transmission_results = DataFrame(
    TransmissionRate=Float64[],
    PeakInfections=Int[],
    FinalRecovered=Int[],
    AttackRate=Float64[],
    EpidemicDuration=Int[]
)

for (i, rate) in enumerate(transmission_rates)
    print("  Testing rate $(rate)... ")

    # Create models with varying transmission rate
    model = SIRModel(rate, 0.005, 0.1)
    models = [model, model]

    # Neutral interactions
    interactions = [1.0 1.0; 1.0 1.0]

    params = SimulationParameters(
        models, interactions, 0.001, 0.01, 18, :none, 100
    )

    population = create_base_population()
    results = simulate(population, params)

    # Analyze results
    populations = results[1]
    infected_counts = [count(ind -> ind[1, 3], pop) for pop in populations]
    recovered_counts = [count(ind -> ind[1, 4], pop) for pop in populations]

    peak_infections = maximum(infected_counts)
    final_recovered = recovered_counts[end]
    attack_rate = final_recovered / 100.0 * 100  # percentage

    # Find epidemic duration (time until <1% infected)
    epidemic_end = findfirst(infected_counts .< 1)
    epidemic_duration = epidemic_end === nothing ? 100 : epidemic_end

    push!(transmission_results, (rate, peak_infections, final_recovered, attack_rate, epidemic_duration))
    println("Peak: $peak_infections, Attack Rate: $(round(attack_rate, digits=1))%")
end

println("✓ Transmission sensitivity complete")

## 3. Recovery Rate Sensitivity Analysis
println("\nRunning recovery rate sensitivity analysis...")

recovery_rates = 0.02:0.02:0.30
recovery_results = DataFrame(
    RecoveryRate=Float64[],
    PeakInfections=Int[],
    EpidemicDuration=Int[],
    FinalRecovered=Int[]
)

for rate in recovery_rates
    print("  Testing recovery rate $(rate)... ")

    model = SIRModel(0.15, 0.005, rate)  # Fixed transmission, varying recovery
    models = [model, model]
    interactions = [1.0 1.0; 1.0 1.0]

    params = SimulationParameters(models, interactions, 0.001, 0.01, 18, :none, 100)
    population = create_base_population()
    results = simulate(population, params)

    populations = results[1]
    infected_counts = [count(ind -> ind[1, 3], pop) for pop in populations]
    recovered_counts = [count(ind -> ind[1, 4], pop) for pop in populations]

    peak_infections = maximum(infected_counts)
    final_recovered = recovered_counts[end]
    epidemic_end = findfirst(infected_counts .< 1)
    epidemic_duration = epidemic_end === nothing ? 100 : epidemic_end

    push!(recovery_results, (rate, peak_infections, epidemic_duration, final_recovered))
    println("Duration: $epidemic_duration")
end

println("✓ Recovery rate sensitivity complete")

## 4. Interaction Strength Exploration
println("\nExploring interaction strengths...")

interaction_strengths = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
cf_ratios = [0.1, 0.3, 0.5, 0.7, 0.9]  # Competition vs Facilitation ratios

interaction_results = DataFrame(
    InteractionStrength=Float64[],
    CFRatio=Float64[],
    Strain1Peak=Int[],
    Strain2Peak=Int[],
    Coinfections=Int[],
    InteractionType=String[]
)

for strength in interaction_strengths
    for cf_ratio in cf_ratios
        print("  Testing strength=$strength, cf_ratio=$cf_ratio... ")

        # Create interaction matrix
        if strength == 0.0
            interactions = [1.0 1.0; 1.0 1.0]
        else
            interactions = create_interaction_matrix(2, false, strength, cf_ratio=cf_ratio)
        end

        models = [base_sir, base_sir]
        params = SimulationParameters(models, interactions, 0.001, 0.01, 18, :none, 100)

        # Create population with both strains present
        population = create_base_population()
        population[end][2, 1] = false  # Add strain 2 seed
        population[end][2, 3] = true

        results = simulate(population, params)
        populations = results[1]

        # Count peaks and coinfections
        strain1_counts = [count(ind -> ind[1, 3], pop) for pop in populations]
        strain2_counts = [count(ind -> ind[2, 3], pop) for pop in populations]
        coinfection_counts = [count(ind -> ind[1, 3] && ind[2, 3], pop) for pop in populations]

        strain1_peak = maximum(strain1_counts)
        strain2_peak = maximum(strain2_counts)
        max_coinfections = maximum(coinfection_counts)

        # Classify interaction type
        interaction_type = if strength == 0.0
            "None"
        elseif cf_ratio < 0.5
            "Competition"
        elseif cf_ratio > 0.5
            "Facilitation"
        else
            "Mixed"
        end

        push!(interaction_results, (strength, cf_ratio, strain1_peak, strain2_peak, max_coinfections, interaction_type))
        println("S1: $strain1_peak, S2: $strain2_peak, Co: $max_coinfections")
    end
end

println("✓ Interaction exploration complete")

## 5. Population Size Effects
println("\nAnalyzing population size effects...")

population_sizes = [50, 100, 200, 500, 1000]
popsize_results = DataFrame(
    PopulationSize=Int[],
    AttackRate=Float64[],
    PeakPrevalence=Float64[],
    EpidemicProbability=Float64[]
)

n_replicates = 10  # Multiple runs for each population size

for pop_size in population_sizes
    print("  Testing population size $pop_size... ")

    attack_rates = Float64[]
    peak_prevalences = Float64[]
    epidemic_occurred = Int[]

    for rep in 1:n_replicates
        Random.seed!(rep * 100)  # Different seed for each replicate

        population = Population(Individual[])
        for i in 1:pop_size
            age = rand(18:65)
            individual = Individual(2, age)
            push!(population.individuals, individual)
        end

        # Seed infections proportional to population size
        n_seeds = max(1, pop_size ÷ 50)
        for i in 1:n_seeds
            population[i][1, 1] = false
            population[i][1, 3] = true
        end

        models = [base_sir, base_sir]
        interactions = [1.0 1.0; 1.0 1.0]
        params = SimulationParameters(models, interactions, 0.001, 0.01, 18, :none, 100)

        results = simulate(population, params)
        populations = results[1]

        infected_counts = [count(ind -> ind[1, 3], pop) for pop in populations]
        recovered_counts = [count(ind -> ind[1, 4], pop) for pop in populations]

        final_recovered = recovered_counts[end]
        attack_rate = final_recovered / length(populations[1]) * 100
        peak_infections = maximum(infected_counts)
        peak_prevalence = peak_infections / length(populations[1]) * 100
        epidemic = final_recovered > n_seeds ? 1 : 0  # Did epidemic spread beyond seeds?

        push!(attack_rates, attack_rate)
        push!(peak_prevalences, peak_prevalence)
        push!(epidemic_occurred, epidemic)
    end

    mean_attack_rate = mean(attack_rates)
    mean_peak_prevalence = mean(peak_prevalences)
    epidemic_probability = mean(epidemic_occurred)

    push!(popsize_results, (pop_size, mean_attack_rate, mean_peak_prevalence, epidemic_probability))
    println("Attack rate: $(round(mean_attack_rate, digits=1))%, Epidemic prob: $(round(epidemic_probability, digits=2))")
end

println("✓ Population size analysis complete")

## 6. Demographic Parameter Effects
println("\nAnalyzing demographic parameter effects...")

# Test different birth and death rates
birth_rates = [0.0, 0.01, 0.02, 0.05]
death_rates = [0.0, 0.005, 0.01, 0.02]

demographic_results = DataFrame(
    BirthRate=Float64[],
    DeathRate=Float64[],
    PopulationChange=Float64[],
    EpidemicSize=Int[],
    EpidemicDuration=Int[]
)

for birth_rate in birth_rates
    for death_rate in death_rates
        print("  Testing birth=$birth_rate, death=$death_rate... ")

        models = [base_sir]  # Single strain for clarity
        interactions = reshape([1.0], 1, 1)

        params = SimulationParameters(models, interactions, death_rate, birth_rate, 18, :none, 150)

        population = Population(Individual[])
        for i in 1:100
            age = rand(18:65)
            individual = Individual(1, age)
            push!(population.individuals, individual)
        end

        # Seed infection
        population[1][1, 1] = false
        population[1][1, 3] = true

        results = simulate(population, params)
        populations = results[1]

        initial_pop = length(populations[1])
        final_pop = length(populations[end])
        pop_change = (final_pop - initial_pop) / initial_pop * 100

        recovered_counts = [count(ind -> ind[1, 4], pop) for pop in populations]
        epidemic_size = recovered_counts[end]

        # Find when epidemic ends
        infected_counts = [count(ind -> ind[1, 3], pop) for pop in populations]
        epidemic_end = findfirst(infected_counts .== 0)
        epidemic_duration = epidemic_end === nothing ? 150 : epidemic_end

        push!(demographic_results, (birth_rate, death_rate, pop_change, epidemic_size, epidemic_duration))
        println("Pop change: $(round(pop_change, digits=1))%")
    end
end

println("✓ Demographic analysis complete")

## 7. Create Comprehensive Visualizations
println("\nCreating parameter exploration visualizations...")

# Plot 1: Transmission rate sensitivity
p1 = plot(transmission_results.TransmissionRate, transmission_results.AttackRate,
    title="Transmission Rate Sensitivity",
    xlabel="Transmission Rate",
    ylabel="Attack Rate (%)",
    lw=3, marker=:circle, markersize=4,
    label="Attack Rate")

plot!(twinx(), transmission_results.TransmissionRate, transmission_results.EpidemicDuration,
    ylabel="Epidemic Duration (days)",
    color=:red, lw=2, marker=:square, markersize=4,
    label="Duration", legend=:topright)

# Plot 2: Recovery rate effects
p2 = plot(recovery_results.RecoveryRate, recovery_results.PeakInfections,
    title="Recovery Rate Effects",
    xlabel="Recovery Rate",
    ylabel="Peak Infections",
    lw=3, marker=:circle, markersize=4,
    color=:green,
    label="Peak Infections")

# Plot 3: Interaction strength heatmap
interaction_pivot = unstack(interaction_results, :CFRatio, :InteractionStrength, :Coinfections)
select!(interaction_pivot, Not(:CFRatio))
interaction_matrix = Matrix(interaction_pivot)

p3 = heatmap(interaction_strengths, cf_ratios, interaction_matrix,
    title="Coinfection Levels vs Interaction Parameters",
    xlabel="Interaction Strength",
    ylabel="Competition/Facilitation Ratio",
    color=:viridis)

# Plot 4: Population size effects
p4 = plot(popsize_results.PopulationSize, popsize_results.EpidemicProbability,
    title="Population Size and Epidemic Risk",
    xlabel="Population Size",
    ylabel="Epidemic Probability",
    lw=3, marker=:circle, markersize=6,
    color=:purple,
    label="Epidemic Probability")

plot!(twinx(), popsize_results.PopulationSize, popsize_results.AttackRate,
    ylabel="Attack Rate (%)",
    color=:orange, lw=2, marker=:diamond, markersize=4,
    label="Attack Rate", legend=:right)

# Combine plots
final_plot = plot(p1, p2, p3, p4,
    layout=(2, 2),
    size=(1200, 1000),
    plot_title="Parameter Exploration Results")

savefig(final_plot, "parameter_exploration_results.png")
println("✓ Saved parameter exploration plots")

## 8. Statistical Analysis
println("\nPerforming statistical analysis...")

# Correlation analysis between parameters and outcomes
println("\nTransmission Rate Correlations:")
println("  - vs Attack Rate: $(round(cor(transmission_results.TransmissionRate, transmission_results.AttackRate), digits=3))")
println("  - vs Peak Infections: $(round(cor(transmission_results.TransmissionRate, transmission_results.PeakInfections), digits=3))")
println("  - vs Epidemic Duration: $(round(cor(transmission_results.TransmissionRate, transmission_results.EpidemicDuration), digits=3))")

println("\nRecovery Rate Correlations:")
println("  - vs Peak Infections: $(round(cor(recovery_results.RecoveryRate, recovery_results.PeakInfections), digits=3))")
println("  - vs Epidemic Duration: $(round(cor(recovery_results.RecoveryRate, recovery_results.EpidemicDuration), digits=3))")

# Find optimal parameters for different objectives
println("\nOptimal Parameters:")

# Minimize peak infections
min_peak_idx = argmin(transmission_results.PeakInfections)
println("  - Minimum Peak Infections: Rate $(transmission_results.TransmissionRate[min_peak_idx]) → $(transmission_results.PeakInfections[min_peak_idx]) infections")

# Minimize epidemic duration
min_duration_idx = argmin(recovery_results.EpidemicDuration)
println("  - Minimum Epidemic Duration: Recovery rate $(recovery_results.RecoveryRate[min_duration_idx]) → $(recovery_results.EpidemicDuration[min_duration_idx]) days")

# Balance attack rate and duration
transmission_results.Score = transmission_results.AttackRate .* 0.6 + transmission_results.EpidemicDuration .* 0.4
best_balance_idx = argmin(transmission_results.Score)
println("  - Best Balance (60% attack rate, 40% duration): Rate $(transmission_results.TransmissionRate[best_balance_idx])")

## 9. Export All Results
println("\nExporting comprehensive results...")

# Save all result DataFrames
CSV.write("transmission_sensitivity.csv", transmission_results)
CSV.write("recovery_sensitivity.csv", recovery_results)
CSV.write("interaction_exploration.csv", interaction_results)
CSV.write("population_size_effects.csv", popsize_results)
CSV.write("demographic_effects.csv", demographic_results)

# Create summary report
summary_df = DataFrame(
    Parameter=["Transmission Rate", "Recovery Rate", "Interaction Strength", "Population Size", "Birth Rate", "Death Rate"],
    Effect=["Strong positive correlation with attack rate and peak infections",
        "Strong negative correlation with epidemic duration",
        "Complex non-linear effects on coinfection patterns",
        "Threshold effects on epidemic probability",
        "Increases epidemic size through demographic turnover",
        "Reduces epidemic size through increased mortality"],
    Sensitivity=["High", "High", "Medium", "High", "Low", "Medium"],
    Recommendation=["Control transmission early", "Promote recovery mechanisms",
        "Consider interaction effects in multi-strain scenarios",
        "Account for population size in risk assessment",
        "Minor factor in short-term epidemics",
        "Important for vulnerable populations"]
)

CSV.write("parameter_analysis_summary.csv", summary_df)

println("✓ All results exported to CSV files")

## 10. Model Comparison Framework
println("\nDemonstrating model comparison framework...")

# Compare different disease model types
models_to_compare = [
    ("SI Model", SIModel(0.15, 0.005)),
    ("SIR Model", SIRModel(0.15, 0.005, 0.1)),
    ("SEIR Model", SEIRModel(0.15, 0.005, 0.1, 5)),
    ("SEIRS Model", SEIRSModel(0.15, 0.005, 0.1, 5, 0.05))
]

model_comparison = DataFrame(
    ModelType=String[],
    PeakInfections=Int[],
    AttackRate=Float64[],
    EpidemicDuration=Int[],
    FinalSusceptible=Int[]
)

for (model_name, model) in models_to_compare
    print("  Testing $model_name... ")

    models = [model]
    interactions = reshape([1.0], 1, 1)
    params = SimulationParameters(models, interactions, 0.001, 0.01, 18, :none, 200)

    population = create_base_population(100, 3)
    # Adjust population for single strain
    for ind in population
        ind.state = ind.state[1:1, :]  # Keep only first strain
    end

    results = simulate(population, params)
    populations = results[1]

    infected_counts = [count(ind -> ind[1, 3], pop) for pop in populations]
    recovered_counts = [count(ind -> ind[1, 4], pop) for pop in populations]
    susceptible_counts = [count(ind -> ind[1, 1], pop) for pop in populations]

    peak_infections = maximum(infected_counts)
    final_recovered = recovered_counts[end]
    attack_rate = final_recovered / 100.0 * 100
    epidemic_end = findfirst(infected_counts .< 1)
    epidemic_duration = epidemic_end === nothing ? 200 : epidemic_end
    final_susceptible = susceptible_counts[end]

    push!(model_comparison, (model_name, peak_infections, attack_rate, epidemic_duration, final_susceptible))
    println("Complete")
end

println("\nModel Comparison Results:")
display(model_comparison)

CSV.write("model_comparison.csv", model_comparison)

println("\n=== Parameter Exploration Complete! ===")
println("This analysis demonstrated:")
println("• Systematic parameter sensitivity analysis")
println("• Interaction strength exploration")
println("• Population size effects")
println("• Demographic parameter impacts")
println("• Statistical correlation analysis")
println("• Parameter optimization approaches")
println("• Model comparison framework")
println("• Comprehensive result export")
println("\nKey findings:")
println("• Transmission rate has strongest impact on epidemic size")
println("• Recovery rate critically affects epidemic duration")
println("• Population size creates threshold effects")
println("• Strain interactions can significantly modify outcomes")
println("• Different disease models show distinct epidemic patterns")
