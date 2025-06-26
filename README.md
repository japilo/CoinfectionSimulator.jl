# CoinfectionSimulator.jl

[![CI](https://github.com/japilo/CoinfectionSimulator.jl/workflows/CI/badge.svg)](https://github.com/japilo/CoinfectionSimulator.jl/actions)
[![codecov](https://codecov.io/gh/japilo/CoinfectionSimulator.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/japilo/CoinfectionSimulator.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://japilo.github.io/CoinfectionSimulator.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://japilo.github.io/CoinfectionSimulator.jl/dev)

A Julia package for simulating multi-strain coinfection dynamics in host populations with imperfect ecological sampling.

The CoinfectionSimulator package includes:
- ✅ Core simulation engine with SI/SIR/SEIR/SEIRS disease models
- ✅ Multi-strain coinfection dynamics with interaction effects  
- ✅ Demographic processes (birth, death, aging)
- ✅ Virtual ecologist sampling with detection errors
- ✅ Interaction matrix generation utilities
- ✅ Comprehensive test suite (47/47 tests passing)
- ✅ Full documentation and examples

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

## Quick Start

```julia
using CoinfectionSimulator
using Random

# Set up initial population (100 individuals, 2 strains)
n_individuals = 100
n_strains = 2
initial_pop = [falses(n_strains, 4) for _ in 1:n_individuals]  # [S,E,I,R] states

# Everyone starts susceptible
for individual in initial_pop
    individual[:, 1] .= true  # Susceptible
end

# Introduce initial infections
initial_pop[1][1, 1] = false  # Individual 1, strain 1: not susceptible
initial_pop[1][1, 3] = true   # Individual 1, strain 1: infected

# Set up simulation parameters
ages = rand(1:50, n_individuals)  # Random ages
interactions = [1.0 0.8; 0.8 1.0]  # Interaction matrix
disease_types = ["si", "si"]  # Both strains use SI model
transmission = [0.3, 0.4]  # Transmission rates

# Run simulation
results = coinfection_simulator(
    initial_pop=initial_pop,
    ages=ages,
    interactions=interactions,
    disease_type=disease_types,
    base_mortality=0.01,
    disease_mortality=[0.02, 0.02],
    fecundity=0.1,
    transmission=transmission,
    time_steps=50,
    age_maturity=18
)

populations, age_vectors = results

# Sample with virtual ecologist
detections = virtual_ecologist_sample(
    virtual_population=populations,
    proportion_sampled=0.3,
    false_positive_rate=0.05,
    false_negative_rate=0.1
)

println("Detected strains over time:")
println(size(detections))  # (time_steps, n_strains)
```

## Main Functions

### `coinfection_simulator`

Simulates multi-strain coinfection dynamics in a host population.

**Key parameters:**
- `initial_pop`: Vector of individual disease state matrices
- `interactions`: Strain interaction matrix
- `disease_type`: Vector of disease model types ("si", "sir", "seir", "seirs")
- `transmission`: Vector of transmission rates
- `time_steps`: Number of simulation time steps

### `virtual_ecologist_sample`

Simulates ecological sampling with imperfect detection.

**Key parameters:**
- `virtual_population`: Simulation results from `coinfection_simulator`
- `proportion_sampled`: Fraction of population sampled
- `false_positive_rate`: Probability of false positive detection
- `false_negative_rate`: Probability of false negative detection

### `prep_interaction_matrix`

Generates interaction matrices for simulation experiments.

**Key parameters:**
- `df`: DataFrame with columns for `interaction_strength`, `cf_ratio`, `priority_effects`, `strains`

## Disease Models

- **SI**: Susceptible → Infected
- **SIR**: Susceptible → Infected → Recovered  
- **SEIR**: Susceptible → Exposed → Infected → Recovered
- **SEIRS**: Susceptible → Exposed → Infected → Recovered → Susceptible

Each strain can use a different disease model within the same simulation.

## Examples

See the `examples/` directory for detailed usage examples:
- Basic single-strain simulation
- Multi-strain coinfection with competition
- Parameter sensitivity analysis
- Virtual ecologist sampling scenarios

## Testing

```julia
using Pkg
Pkg.test("CoinfectionSimulator")
```

## Citation

If you use this package in your research, please cite:

```
Pilowsky, J. (2025). CoinfectionSimulator.jl: A Julia package for multi-strain coinfection dynamics.
```

## License

MIT License - see LICENSE file for details.
