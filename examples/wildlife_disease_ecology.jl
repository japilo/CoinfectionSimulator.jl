# Wildlife Disease Ecology Example: Multiple Infections in Animal Populations
#
# This example demonstrates CoinfectionSimulator.jl for wildlife disease ecology
# research with directly transmitted parasites and pathogens, including:
# - Parasites and pathogens with density-dependent transmission
# - Animal population dynamics with birth, death, and aging
# - Strain interactions in coinfected individuals
# - Wildlife sampling challenges and detection limitations

using CoinfectionSimulator
using Random
using Plots
using DataFrames
using Statistics
using CSV
Random.seed!(2025)

## 1. Define Parasite Models for Wildlife Host
# In this example, we imagine a host population of wildlife that can be infected with four different parasites/pathogens.
# We will consider the time step of the simulation to mean a month in the real world (you can interpret the time steps of the simulation as you like, though a year is likely too long.)

# Strain 1: Acute viral infection
# High transmission, short infectious period, immunity develops
acute_virus = SIRModel(
    0.08,    # Transmission rate per contact per month.
    0.012,   # Disease-induced mortality
    0.25     # Recovery rate (4 month average infectious period)
)

# Strain 2: Chronic bacterial infection
# Moderate transmission, animals remain infectious, no recovery
chronic_bacteria = SIModel(
    0.05,    # Lower transmission rate
    0.008    # Lower mortality but persistent
)

# Strain 3: Parasitic worm with prepatent period
# Transmission through contact, latent period before becoming infectious
intestinal_parasite = SEIRModel(
    0.06,    # Moderate transmission
    0.005,   # Low mortality
    0.15,    # Recovery/clearance rate
    6        # Prepatent period (6 months before infectious)
)

# Strain 4: Ectoparasite
ectoparasite = SEIRSModel(
    0.10,    # High transmission
    0.003,   # Very low mortality
    0.20,    # Recovery rate
    2,       # Short latency
    0.08     # Loss of immunity/resistance
)

disease_models = [acute_virus, chronic_bacteria, intestinal_parasite, ectoparasite]

## 2. Define Strain Interactions Based on Immunological and Ecological Mechanisms

# Interactions based on biological mechanisms:
# - Cross-reactive immunity between related strains
# - Immunosuppression by chronic infections
# - Competition for host resources
# - Facilitation through immune system damage
# - Priority effects (matrix is asymmetric, it matters which strain infects the host first)

interactions = [
    1.0 0.7 1.2 0.9;    # Acute virus: competes with bacteria, facilitates internal parasite
    1.3 1.0 1.1 1.0;    # Chronic bacteria: facilitates others (immunosuppression)
    0.8 0.9 1.0 1.1;    # Intestinal parasite: some competition, slight facilitation of ectoparasite
    1.1 1.0 0.9 1.0     # Ectoparasite: facilitates virus, competes with intestinal parasite
]

## 3. Initialize Wildlife Population with Realistic Demographics

function create_wildlife_population(n_animals=300)
    population = Population(Individual[])

    # Age structure typical of wildlife populations
    # Most animals are young, fewer survive to old age
    for i in 1:n_animals
        # Age in time steps (represents months)
        age = if rand() < 0.45
            rand(0:6)       # 45% juveniles (0-6 months)
        elseif rand() < 0.35
            rand(7:18)      # 35% young adults (7-18 months)
        elseif rand() < 0.15
            rand(19:36)     # 15% mature adults (19-36 months)
        else
            rand(37:60)     # 5% old animals (3-5 years)
        end

        individual = Individual(4, age)  # 4 strains
        push!(population.individuals, individual)
    end

    return population
end

wildlife_population = create_wildlife_population(300)

# Create a new population with all individuals initially susceptible
wildlife_population = create_wildlife_population(300)

# Now let's introduce infections in a safe way that ensures each strain
# is in exactly one state per individual

# Chronic bacteria most common (persists in population)
chronic_initial = 8
for i in 1:chronic_initial
    # Make infected with chronic bacteria (strain 2)
    # First set susceptible to false
    wildlife_population.individuals[i].state[2, 1] = false
    # Then set infected to true
    wildlife_population.individuals[i].state[2, 3] = true
end

# Acute virus (recent outbreak)
virus_initial = 3
for i in (chronic_initial+1):(chronic_initial+virus_initial)
    # Make infected with acute virus (strain 1)
    wildlife_population.individuals[i].state[1, 1] = false
    wildlife_population.individuals[i].state[1, 3] = true
end

# Intestinal parasites (some latent infections from previous exposure)
parasite_initial = 5
for i in (chronic_initial+virus_initial+1):(chronic_initial+virus_initial+parasite_initial)
    # Make exposed to intestinal parasite (strain 3)
    wildlife_population.individuals[i].state[3, 1] = false
    wildlife_population.individuals[i].state[3, 2] = true  # Exposed (prepatent)
end

# Ectoparasites (few, transmission just starting)
ecto_initial = 2
for i in (chronic_initial+virus_initial+parasite_initial+1):(chronic_initial+virus_initial+parasite_initial+ecto_initial)
    # Make infected with ectoparasite (strain 4)
    wildlife_population.individuals[i].state[4, 1] = false
    wildlife_population.individuals[i].state[4, 3] = true
end

## 4. Run Multi-Season Simulation

# Simulate 2 years (24 time steps = months)
params = SimulationParameters(
    disease_models,
    interactions,
    0.008,   # Base mortality (predation, accidents, etc.)
    0.05,    # Fecundity (birth rate)
    6,       # Reproductive maturity at 6 months
    :random, # Random introductions (migration, spillover events)
    24       # 2 years of monthly time steps
)

populations = simulate(wildlife_population, params)

## 5. Analyze Disease Dynamics and Population Health

n_timesteps = length(populations)
n_strains = 4
strain_names = ["Acute Virus", "Chronic Bacteria", "Intestinal Parasite", "Ectoparasite"]

# Track all disease states
susceptible = zeros(Int, n_timesteps, n_strains)
exposed = zeros(Int, n_timesteps, n_strains)
infected = zeros(Int, n_timesteps, n_strains)
recovered = zeros(Int, n_timesteps, n_strains)
total_animals = zeros(Int, n_timesteps)

# Track population health metrics
total_infected_any = zeros(Int, n_timesteps)
multiple_infections = zeros(Int, n_timesteps)
uninfected_animals = zeros(Int, n_timesteps)

for t in 1:n_timesteps
    pop = populations[t]
    total_animals[t] = length(pop.individuals)

    for i in 1:length(pop.individuals)
        animal = pop.individuals[i]
        # Count disease states for each strain
        for strain in 1:n_strains
            if animal.state[strain, 1]
                susceptible[t, strain] += 1
            elseif animal.state[strain, 2]
                exposed[t, strain] += 1
            elseif animal.state[strain, 3]
                infected[t, strain] += 1
            elseif animal.state[strain, 4]
                recovered[t, strain] += 1
            end
        end

        # Population health metrics
        infections_count = sum(animal.state[strain, 3] for strain in 1:n_strains)
        any_infection = any(animal.state[strain, 3] for strain in 1:n_strains)

        if any_infection
            total_infected_any[t] += 1
        else
            uninfected_animals[t] += 1
        end

        if infections_count >= 2
            multiple_infections[t] += 1
        end
    end
end

# Calculate key epidemiological metrics

for strain in 1:n_strains
    peak_infected = maximum(infected[:, strain])
    peak_month = findfirst(infected[:, strain] .== peak_infected)
    total_ever_infected = maximum(infected[:, strain] .+ recovered[:, strain])
    prevalence_final = infected[end, strain] / total_animals[end] * 100
end

# Coinfection analysis
max_multiple = maximum(multiple_infections)
avg_multiple = mean(multiple_infections)

# Population health
final_uninfected = uninfected_animals[end] / total_animals[end] * 100

## 6. Age-structured Disease Analysis

# Examine final population by age groups
final_pop = populations[end]
age_groups = [(0, 6, "Juvenile"), (7, 18, "Young Adult"), (19, 36, "Mature Adult"), (37, 100, "Old")]

age_disease_df = DataFrame(
    AgeGroup=String[],
    Count=Int[],
    Uninfected=Int[],
    SingleInfection=Int[],
    MultipleInfections=Int[]
)

for (min_age, max_age, group_name) in age_groups
    group_animals = [animal for animal in final_pop.individuals if min_age <= animal.age <= max_age]
    count = length(group_animals)

    if count > 0
        uninfected = sum(!any(animal.state[strain, 3] for strain in 1:4) for animal in group_animals)
        single_inf = sum(sum(animal.state[strain, 3] for strain in 1:4) == 1 for animal in group_animals)
        multiple_inf = sum(sum(animal.state[strain, 3] for strain in 1:4) >= 2 for animal in group_animals)

        push!(age_disease_df, (group_name, count, uninfected, single_inf, multiple_inf))
    end
end

display(age_disease_df)

## 7. Wildlife Sampling Simulation
# Here we will explore different possible scenarios of disease surveillance in wildlife and how they affect detection of strains.

# Simulate different sampling strategies used in wildlife research
sampling_strategies = [
    ("Capture-Mark-Recapture", 0.15, 0.02, 0.08),  # 15% sampled, low errors (blood samples)
    ("Hunter Harvest", 0.08, 0.05, 0.12),          # 8% sampled, moderate errors (field necropsies)
    ("Mortality Surveillance", 0.03, 0.10, 0.20)    # 3% sampled, high errors (found dead animals)
]

surveillance_results = DataFrame(
    Strategy=String[],
    Strain=String[],
    TruePositiveMonths=Int[],
    DetectedMonths=Int[],
    SensitivityPercent=Float64[]
)

for (strategy_name, prop_sampled, fp_rate, fn_rate) in sampling_strategies

    sampling_params = SamplingParameters(prop_sampled, fp_rate, fn_rate)
    detection_matrix = sample_populations(populations, sampling_params)

    for strain in 1:n_strains
        # Months where strain was actually present
        true_positive_months = sum((infected[:, strain] ./ total_animals) .> 0)

        # Months where strain was detected by sampling
        detected_months = sum(detection_matrix[:, strain])

        sensitivity = detected_months / max(1, true_positive_months) * 100

        push!(surveillance_results, (strategy_name, strain_names[strain],
            true_positive_months, detected_months, sensitivity))
    end
end

display(surveillance_results)

## 8. Create Wildlife Disease Visualizations

months = 1:n_timesteps

# Plot 1: Multi-strain prevalence over time
p1 = plot(title="Wildlife Disease Dynamics Over 2 Years",
    xlabel="Month", ylabel="Infected Animals")
colors = [:red, :blue, :green, :orange]
for strain in 1:n_strains
    plot!(p1, months, infected[:, strain],
        label=strain_names[strain], lw=2, color=colors[strain])
end

# Plot 2: Population health metrics
p2 = plot(title="Population Health Status", xlabel="Month", ylabel="Number of Animals")
plot!(p2, months, uninfected_animals, label="Uninfected", lw=2, color=:green)
plot!(p2, months, total_infected_any, label="Any Infection", lw=2, color=:red)
plot!(p2, months, multiple_infections, label="Multiple Infections", lw=2, color=:purple)

# Plot 3: Age-disease relationship
age_groups_plot = age_disease_df.AgeGroup
infected_by_age = age_disease_df.SingleInfection .+ age_disease_df.MultipleInfections
p3 = bar(age_groups_plot, infected_by_age,
    title="Disease Burden by Age Group", ylabel="Infected Animals",
    color=[:lightblue, :orange, :red, :darkred], legend=false)

# Combine plots
final_plot = plot(p1, p2, p3, layout=(3, 1), size=(1400, 1000),
    plot_title="Wildlife Disease Ecology: Multi-strain Dynamics")

savefig(final_plot, "wildlife_disease_ecology.png")

## 9. Export Wildlife Disease Data

# Disease time series
disease_timeseries = DataFrame(
    Month=repeat(1:n_timesteps, 4),
    Strain=repeat(strain_names, inner=n_timesteps),
    Susceptible=vcat([susceptible[:, i] for i in 1:4]...),
    Exposed=vcat([exposed[:, i] for i in 1:4]...),
    Infected=vcat([infected[:, i] for i in 1:4]...),
    Recovered=vcat([recovered[:, i] for i in 1:4]...),
    TotalAnimals=repeat(total_animals, 4)
)

CSV.write("wildlife_disease_timeseries.csv", disease_timeseries)

# Surveillance results
CSV.write("wildlife_surveillance_results.csv", surveillance_results)

# Age-disease patterns
CSV.write("wildlife_age_disease_patterns.csv", age_disease_df)
