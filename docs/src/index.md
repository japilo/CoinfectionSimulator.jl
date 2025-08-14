```@meta
CurrentModule = CoinfectionSimulator
```

# CoinfectionSimulator.jl

Documentation for [CoinfectionSimulator.jl](https://github.com/japilo/CoinfectionSimulator.jl).

A Julia package for simulating multi-strain coinfection dynamics in host populations with imperfect ecological sampling.

There are many ways to simulate compartmental disease models in Julia. What CoinfectionSimulator.jl brings to the table is community ecology interactions between pathogens that circulate in a host population. The simulator works off of the assumption that pathogens can affect the probability of other pathogens successfully infecting the host, either by facilitating or inhibiting their establishment. The user can specify all possible pairwise interactions between pathogens using an interaction matrix, which can include or exclude priority effects. In its current form, the simulator assumes density-dependent transmission.

## Features

- **Multi-strain epidemiological models**: Supports SI, SIR, SEIR, and SEIRS disease models
- **Coinfection dynamics**: Models strain interactions (competition, facilitation, priority effects)
- **Population dynamics**: Includes birth, death, and aging processes
- **Virtual ecologist sampling**: Simulates imperfect ecological surveys with false positives/negatives
- **Flexible parameterization**: Easy setup of interaction matrices and simulation parameters

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/japilo/CoinfectionSimulator.jl")
```

## Performance

Here are some performance benchmarks on a 2024 MacBook Pro with an M4 chip and 10 cores, using reasonable simulation parameters (similar base mortality and fecundity in the host population, low virulence strains):

### Population Size (100 timesteps)

| Population Size | Time     | Memory   |
|-----------------|----------|----------|
| 100             | 0.85 ms  | 5.2 MB   |
| 250             | 2.27 ms  | 13.5 MB  |
| 500             | 5.15 ms  | 26.7 MB  |
| 1000            | 18.83 ms | 54.4 MB  |
| 2000            | 41.06 ms | 112.0 MB |

### Time Steps (500 individuals)

| Time Steps | Time      | Memory    |
|------------|-----------|-----------|
| 50         | 2.32 ms   | 11.6 MB   |
| 100        | 4.60 ms   | 23.6 MB   |
| 200        | 10.07 ms  | 45.5 MB   |
| 300        | 23.48 ms  | 64.9 MB   |
| 400        | 43.24 ms  | 100.2 MB  |
| 500        | 46.70 ms  | 122.2 MB  |

## Quick Start

```julia
using CoinfectionSimulator
using Random

# Set up initial population (100 individuals, 2 strains)
n_individuals = 100
n_strains = 2

# Create disease models
models = [SIModel(0.3, 0.02), SIModel(0.4, 0.02)]

# Create initial population
pop = Population([Individual(BitMatrix([true false true false]), 20) for _ in 1:100])

# Set up simulation parameters
params = SimulationParameters(
    models,                  # Disease models
    [1.0 0.8; 0.8 1.0],      # Interaction matrix
    0.01,                    # Base mortality
    0.1,                     # Fecundity
    30,                      # Age of maturity
    :simultaneous,           # Introduction type
    50,                      # Time steps
    :frequency               # Transmission type (:frequency or :density)
)

# Run simulation
populations, age_vectors = simulate(population, params)

# Sample with virtual ecologist
sampling_params = SamplingParameters(0.3, 0.05, 0.1)
detections = sample_populations(populations, sampling_params)

println("Detected strains over time:")
println(size(detections))  # (time_steps, n_strains)
```

## Main Types and Functions

### Types

#### Disease Models
- `SIModel`: Susceptible → Infected
- `SIRModel`: Susceptible → Infected → Recovered
- `SEIRModel`: Susceptible → Exposed → Infected → Recovered
- `SEIRSModel`: Susceptible → Exposed → Infected → Recovered → Susceptible

#### Core Types
- `Individual`: Represents a single host with an age in days and disease states for multiple strains
- `Population`: Collection of individuals
- `SimulationParameters`: Configuration for a simulation run
- `SamplingParameters`: Configuration for virtual ecologist sampling

### Primary Functions

#### `simulate`

Simulates multi-strain coinfection dynamics in a host population.

**Parameters:**
- `initial_population`: Population of individuals
- `params`: SimulationParameters object

#### `sample_populations`

Simulates ecological sampling with imperfect detection.

**Parameters:**
- `populations`: Vector of populations across time steps
- `params`: SamplingParameters object

#### `create_interaction_matrix`

Generates interaction matrices for simulation experiments.

**Parameters:**
- `df`: DataFrame with columns for `interaction_strength`, `cf_ratio`, `priority_effects`, `strains`
OR
- `strains`: Number of strains
- `priority_effects`: Whether to use asymmetric interactions
- `interaction_strength`: Strength of strain interactions
- `cf_ratio`: Ratio of competition to facilitation

## Disease Models

- **SI**: Susceptible → Infected
- **SIR**: Susceptible → Infected → Recovered
- **SEIR**: Susceptible → Exposed → Infected → Recovered
- **SEIRS**: Susceptible → Exposed → Infected → Recovered → Susceptible

Each strain can use a different disease model within the same simulation.

```@index
```
