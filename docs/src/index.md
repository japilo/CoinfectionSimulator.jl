```@meta
CurrentModule = CoinfectionSimulator
```

# CoinfectionSimulator.jl

Documentation for [CoinfectionSimulator.jl](https://github.com/japilo/CoinfectionSimulator.jl).

## Overview

CoinfectionSimulator.jl is a Julia package for simulating multi-strain coinfection dynamics in host populations with imperfect ecological sampling.

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

# Create interaction matrix
interactions = prep_interaction_matrix(3, "symmetric", 0.2)

# Set up initial population
initial_pop = [Bool[1 0 0 0; 1 0 0 0; 1 0 0 0] for _ in 1:100]

# Run simulation
result = coinfection_simulator(
    initial_pop=initial_pop,
    ages=rand(1:50, 100),
    interactions=interactions,
    disease_type=["SEIR", "SEIR", "SEIR"],
    base_mortality=0.01,
    disease_mortality=[0.02, 0.02, 0.02],
    fecundity=0.1,
    transmission=[0.3, 0.3, 0.3],
    time_steps=50,
    age_maturity=18,
    latency=[2, 2, 2],
    recovery=[0.1, 0.1, 0.1]
)
```

```@index
```
