```@meta
CurrentModule = CoinfectionSimulator
```

# API Reference

## Disease Model Types

```@docs
DiseaseModel
SIModel
SIRModel
SEIRModel
SEIRSModel
```

## Population Types

```@docs
Individual
Population
```

## Parameter Types

```@docs
SimulationParameters
SamplingParameters
```

## Main Functions

```@docs
simulate
sample_populations
create_interaction_matrix
```

## Helper Functions

```@docs
age_population!
apply_base_mortality!
breeding!
copy_individual
copy_population
introduce_infections!
process_disease_mortality!
process_disease_dynamics!
process_strain!
process_exposures!
process_infections!
process_recovery!
process_latent_infections!
process_immunity_loss!
remove_dead_individuals!
```

## Index

```@index
```
