# Examples

## Basic Single-Strain Simulation

```julia
using CoinfectionSimulator
using Random

Random.seed!(42)

# Single strain, SI model
initial_pop = [Bool[1 0 0 0] for _ in 1:50]  # 50 individuals, all susceptible
initial_pop[1] = Bool[0 0 1 0]  # First individual infected

result = coinfection_simulator(
    initial_pop=initial_pop,
    ages=rand(1:30, 50),
    interactions=reshape([1.0], 1, 1),
    disease_type=["SI"],
    base_mortality=0.01,
    disease_mortality=[0.05],
    fecundity=0.1,
    transmission=[0.3],
    time_steps=20,
    age_maturity=5
)

populations, ages = result
println("Final population size: ", length(populations[end]))
```

## Multi-Strain Coinfection

```julia
using CoinfectionSimulator

# Two competing strains
n_individuals = 100
initial_pop = [Bool[1 0 0 0; 1 0 0 0] for _ in 1:n_individuals]

# Introduce both strains
initial_pop[1] = Bool[0 0 1 0; 1 0 0 0]  # Individual 1: strain 1 infected
initial_pop[2] = Bool[1 0 0 0; 0 0 1 0]  # Individual 2: strain 2 infected

# Competitive interaction matrix
interactions = [1.0 0.7; 0.7 1.0]  # Mutual inhibition

result = coinfection_simulator(
    initial_pop=initial_pop,
    ages=rand(1:50, n_individuals),
    interactions=interactions,
    disease_type=["SIR", "SIR"],
    base_mortality=0.01,
    disease_mortality=[0.03, 0.03],
    fecundity=0.15,
    transmission=[0.4, 0.4],
    time_steps=100,
    age_maturity=10,
    recovery=[0.1, 0.1]
)
```

## Virtual Ecologist Sampling

```julia
# Using results from previous simulation
sampled_data = virtual_ecologist_sample(
    virtual_population=result[1],
    proportion_sampled=0.3,
    false_positive_rate=0.05,
    false_negative_rate=0.15
)

println("Detection matrix dimensions: ", size(sampled_data))
println("Total detections: ", sum(sampled_data))
```

## Parameter Study Setup

```julia
using DataFrames

# Create parameter space
df = DataFrame(
    interaction_strength = [0.1, 0.2, 0.3],
    cf_ratio = [0.8, 1.0, 1.2],
    priority_effects = [true, false, true],
    strains = [2, 3, 2]
)

# Generate interaction matrices
matrices = prep_interaction_matrix(df=df)

for (i, matrix) in enumerate(matrices)
    println("Matrix $i:")
    display(matrix)
    println()
end
```
