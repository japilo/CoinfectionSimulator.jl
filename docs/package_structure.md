# Documentation

## Package Structure

```
CoinfectionSimulator/
├── Project.toml              # Package metadata and dependencies
├── README.md                 # Main documentation
├── src/                      # Source code
│   ├── CoinfectionSimulator.jl  # Main package module
│   ├── simulator.jl          # Core simulation functions
│   ├── sampling.jl           # Virtual ecologist sampling
│   ├── data_prep.jl          # Data preparation utilities
│   └── utils.jl              # Disease model utilities
├── test/                     # Test suite
│   ├── runtests.jl           # Main test runner
│   ├── test_simulator.jl     # Simulator tests
│   ├── test_sampling.jl      # Sampling tests
│   └── test_data_prep.jl     # Data preparation tests
├── examples/                 # Example scripts
│   └── basic_example.jl      # Basic usage example
└── docs/                     # Documentation directory
    └── package_structure.md  # This file
```

## Core Functions

### `coinfection_simulator`
**Location**: `src/simulator.jl`

Main simulation function that models multi-strain coinfection dynamics in a host population.

**Key Features**:
- Multiple disease compartment models (SI, SIR, SEIR, SEIRS)
- Strain interactions (competition, facilitation, priority effects)
- Host population dynamics (birth, death, aging)
- Flexible strain introduction timing

### `virtual_ecologist_sample`
**Location**: `src/sampling.jl`

Simulates ecological field sampling with imperfect detection capabilities.

**Key Features**:
- Configurable sampling proportion
- False positive and false negative error rates
- Strain-level detection outcomes

### `prep_interaction_matrix`
**Location**: `src/data_prep.jl`

Generates interaction matrices for simulation experiments.

**Key Features**:
- Symmetric or asymmetric matrices
- Parameterized interaction strength
- Competitive vs facilitative interactions

## Disease Compartment Models

### SI Model
- **States**: Susceptible (S) → Infected (I)
- **Parameters**: transmission rate
- **Use case**: Chronic infections without recovery

### SIR Model  
- **States**: Susceptible (S) → Infected (I) → Recovered (R)
- **Parameters**: transmission rate, recovery rate
- **Use case**: Acute infections with permanent immunity

### SEIR Model
- **States**: Susceptible (S) → Exposed (E) → Infected (I) → Recovered (R)
- **Parameters**: transmission rate, latency period, recovery rate
- **Use case**: Infections with incubation period

### SEIRS Model
- **States**: S → E → I → R → S (waning immunity)
- **Parameters**: transmission rate, latency, recovery rate, immunity loss rate
- **Use case**: Infections with temporary immunity

## Data Structures

### Individual Representation
Each individual is represented by an n×4 Boolean matrix where:
- **n**: number of strains
- **Columns**: [Susceptible, Exposed, Infected, Recovered]
- **Rows**: one per strain

### Population Format
- **Type**: `Vector{Matrix{Bool}}`
- **Structure**: Vector of individual matrices
- **Time series**: `Vector{Vector{Matrix{Bool}}}` (time → individuals → state matrix)

### Interaction Matrix
- **Type**: `Matrix{Float64}`
- **Dimensions**: n×n (strains × strains)
- **Diagonal**: Always 1.0 (self-interaction)
- **Off-diagonal**: Interaction coefficients (>1: facilitation, <1: competition)

## Testing

Run the full test suite:
```julia
using Pkg
Pkg.test("CoinfectionSimulator")
```

Run specific test modules:
```julia
# In the package directory
julia --project=test/test_simulator.jl
```

## Examples

See `examples/basic_example.jl` for a complete working example demonstrating:
- Population setup
- Multi-strain simulation
- Visualization

## Dependencies

- **DataFrames.jl**: Data manipulation
- **Distributions.jl**: Statistical distributions
- **StatsBase.jl**: Statistical functions and sampling
- **LinearAlgebra.jl**: Matrix operations
- **Random.jl**: Random number generation

## Performance Notes

- Large populations (>10,000 individuals) may require significant memory
- Long simulations (>1,000 time steps) can be computationally intensive
- Consider reducing strain numbers or using smaller test populations for development
- Simulation is not currently parallelized by default but could benefit from multi-threading

## Extending the Package

### Performance Optimization
- Consider using `StaticArrays.jl` for individual matrices
- Implement multi-threading for population processing
- Add memory-efficient streaming for large time series
