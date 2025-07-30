"""
    DiseaseModel

Abstract supertype for all disease models in the simulation.
"""
abstract type DiseaseModel end

"""
    SIModel <: DiseaseModel

Susceptible-Infected (SI) disease model.

# Fields
- `transmission::Float64`: Transmission probability
- `mortality::Float64`: Additional mortality rate due to infection
"""
struct SIModel <: DiseaseModel
    transmission::Float64
    mortality::Float64

    function SIModel(transmission::Float64, mortality::Float64)
        0 ≤ transmission ≤ 1 || throw(ArgumentError("Transmission rate must be between 0 and 1"))
        0 ≤ mortality ≤ 1 || throw(ArgumentError("Mortality rate must be between 0 and 1"))
        new(transmission, mortality)
    end
end

"""
    SIRModel <: DiseaseModel

Susceptible-Infected-Recovered (SIR) disease model.

# Fields
- `transmission::Float64`: Transmission probability
- `mortality::Float64`: Additional mortality rate due to infection
- `recovery::Float64`: Recovery probability
"""
struct SIRModel <: DiseaseModel
    transmission::Float64
    mortality::Float64
    recovery::Float64

    function SIRModel(transmission::Float64, mortality::Float64, recovery::Float64)
        0 ≤ transmission ≤ 1 || throw(ArgumentError("Transmission rate must be between 0 and 1"))
        0 ≤ mortality ≤ 1 || throw(ArgumentError("Mortality rate must be between 0 and 1"))
        0 ≤ recovery ≤ 1 || throw(ArgumentError("Recovery rate must be between 0 and 1"))
        new(transmission, mortality, recovery)
    end
end

"""
    SEIRModel <: DiseaseModel

Susceptible-Exposed-Infected-Recovered (SEIR) disease model.

# Fields
- `transmission::Float64`: Transmission probability
- `mortality::Float64`: Additional mortality rate due to infection
- `recovery::Float64`: Recovery probability
- `latency::Int`: Number of time steps in exposed state
"""
struct SEIRModel <: DiseaseModel
    transmission::Float64
    mortality::Float64
    recovery::Float64
    latency::Int

    function SEIRModel(transmission::Float64, mortality::Float64, recovery::Float64, latency::Int)
        0 ≤ transmission ≤ 1 || throw(ArgumentError("Transmission rate must be between 0 and 1"))
        0 ≤ mortality ≤ 1 || throw(ArgumentError("Mortality rate must be between 0 and 1"))
        0 ≤ recovery ≤ 1 || throw(ArgumentError("Recovery rate must be between 0 and 1"))
        latency > 0 || throw(ArgumentError("Latency period must be positive"))
        new(transmission, mortality, recovery, latency)
    end
end

"""
    SEIRSModel <: DiseaseModel

Susceptible-Exposed-Infected-Recovered-Susceptible (SEIRS) disease model.

# Fields
- `transmission::Float64`: Transmission probability
- `mortality::Float64`: Additional mortality rate due to infection
- `recovery::Float64`: Recovery probability
- `latency::Int`: Number of time steps in exposed state
- `immunity_loss::Float64`: Probability of losing immunity
"""
struct SEIRSModel <: DiseaseModel
    transmission::Float64
    mortality::Float64
    recovery::Float64
    latency::Int
    immunity_loss::Float64

    function SEIRSModel(transmission::Float64, mortality::Float64, recovery::Float64,
        latency::Int, immunity_loss::Float64)
        0 ≤ transmission ≤ 1 || throw(ArgumentError("Transmission rate must be between 0 and 1"))
        0 ≤ mortality ≤ 1 || throw(ArgumentError("Mortality rate must be between 0 and 1"))
        0 ≤ recovery ≤ 1 || throw(ArgumentError("Recovery rate must be between 0 and 1"))
        latency > 0 || throw(ArgumentError("Latency period must be positive"))
        0 ≤ immunity_loss ≤ 1 || throw(ArgumentError("Immunity loss rate must be between 0 and 1"))
        new(transmission, mortality, recovery, latency, immunity_loss)
    end
end

"""
    Individual

Represents a single individual in the population.

# Fields
- `state::Matrix{Bool}`: Disease states matrix with rows for each strain and columns for [S,E,I,R]
- `age::Int`: Age of the individual
"""
mutable struct Individual
    state::BitMatrix
    age::Int

    function Individual(state::AbstractMatrix{Bool}, age::Int)
        size(state, 2) == 4 || throw(ArgumentError("State matrix must have 4 columns for [S,E,I,R]"))
        age >= 0 || throw(ArgumentError("Age must be non-negative"))

        # Validate that each strain is in exactly one state
        for i in 1:size(state, 1)
            sum(state[i, :]) == 1 || throw(ArgumentError("Each strain must be in exactly one state"))
        end

        new(BitMatrix(state), age)
    end
end

# Constructors for creating new susceptible individuals
Individual(n_strains::Int, age::Int) = Individual(hcat(trues(n_strains), falses(n_strains, 3)), age)

# Make individuals indexable directly to their state
Base.getindex(ind::Individual, args...) = getindex(ind.state, args...)
Base.setindex!(ind::Individual, v, args...) = setindex!(ind.state, v, args...)
Base.size(ind::Individual, args...) = size(ind.state, args...)
Base.copy(ind::Individual) = Individual(copy(ind.state), ind.age)
Base.deepcopy(ind::Individual) = Individual(deepcopy(ind.state), ind.age)

"""
    Population

Represents the entire population in the simulation.

# Fields
- `individuals::Vector{Individual}`: Collection of individuals
"""
struct Population
    individuals::Vector{Individual}

    function Population(individuals::Vector{Individual})
        if !isempty(individuals)
            n_strains = size(individuals[1], 1)
            all(ind -> size(ind, 1) == n_strains, individuals) ||
                throw(ArgumentError("All individuals must have the same number of strains"))
        end
        new(individuals)
    end
end

# Constructor to convert legacy format
function Population(matrices::Vector{<:AbstractMatrix{Bool}}, ages::Vector{Int})
    length(matrices) == length(ages) ||
        throw(ArgumentError("Number of matrices must match number of ages"))
    Population([Individual(m, a) for (m, a) in zip(matrices, ages)])
end

# Make population indexable and iterable
Base.getindex(pop::Population, i) = pop.individuals[i]
Base.setindex!(pop::Population, v, i) = (pop.individuals[i] = v)
Base.iterate(pop::Population, state...) = iterate(pop.individuals, state...)
Base.length(pop::Population) = length(pop.individuals)
Base.size(pop::Population) = (length(pop.individuals),)
Base.isempty(pop::Population) = isempty(pop.individuals)
Base.copy(pop::Population) = Population([copy(ind) for ind in pop.individuals])
Base.deepcopy(pop::Population) = Population([deepcopy(ind) for ind in pop.individuals])

# Additional methods for population operations
function n_strains(pop::Population)
    isempty(pop) ? 0 : size(pop[1], 1)
end

function add_individual!(pop::Population, ind::Individual)
    if !isempty(pop)
        size(ind, 1) == n_strains(pop) ||
            throw(ArgumentError("New individual must have the same number of strains"))
    end
    push!(pop.individuals, ind)
end

"""
    SimulationParameters

Parameters for a coinfection simulation.

# Fields
- `models::Vector{<:DiseaseModel}`: Disease models for each strain
- `interactions::Matrix{Float64}`: Strain interaction matrix
- `base_mortality::Float64`: Background mortality rate
- `fecundity::Float64`: Mean number of offspring per mature individual
- `age_maturity::Int`: Age at which individuals can reproduce
- `introduction::Symbol`: How strains are introduced (:simultaneous, :random, or :none)
- `time_steps::Int`: Number of time steps to simulate
"""
struct SimulationParameters
    models::Vector{<:DiseaseModel}
    interactions::Matrix{Float64}
    base_mortality::Float64
    fecundity::Float64
    age_maturity::Int
    introduction::Symbol
    time_steps::Int

    function SimulationParameters(
        models::Vector{<:DiseaseModel},
        interactions::Matrix{Float64},
        base_mortality::Float64,
        fecundity::Float64,
        age_maturity::Int,
        introduction::Symbol,
        time_steps::Int
    )
        n_strains = length(models)
        size(interactions) == (n_strains, n_strains) ||
            throw(ArgumentError("Interaction matrix size must match number of disease models"))

        0 ≤ base_mortality ≤ 1 || throw(ArgumentError("Base mortality must be between 0 and 1"))
        fecundity ≥ 0 || throw(ArgumentError("Fecundity must be non-negative"))
        age_maturity > 0 || throw(ArgumentError("Age of maturity must be positive"))
        introduction in (:simultaneous, :random, :none) ||
            throw(ArgumentError("Introduction must be :simultaneous, :random, or :none"))
        time_steps ≥ 1 || throw(ArgumentError("Time steps must be positive"))

        # Check that base_mortality + disease_mortality <= 1
        for model in models
            base_mortality + model.mortality ≤ 1 ||
                throw(ArgumentError("Base mortality + disease mortality must not exceed 1"))
        end

        new(models, interactions, base_mortality, fecundity, age_maturity, introduction, time_steps)
    end
end

"""
    SamplingParameters

Parameters for virtual ecologist sampling.

# Fields
- `proportion_sampled::Float64`: Proportion of the population sampled
- `false_positive_rate::Float64`: Probability of false positive detection
- `false_negative_rate::Float64`: Probability of false negative detection
"""
struct SamplingParameters
    proportion_sampled::Float64
    false_positive_rate::Float64
    false_negative_rate::Float64

    function SamplingParameters(
        proportion_sampled::Float64,
        false_positive_rate::Float64,
        false_negative_rate::Float64
    )
        0 ≤ proportion_sampled ≤ 1 || throw(ArgumentError("Proportion sampled must be between 0 and 1"))
        0 ≤ false_positive_rate ≤ 1 || throw(ArgumentError("False positive rate must be between 0 and 1"))
        0 ≤ false_negative_rate ≤ 1 || throw(ArgumentError("False negative rate must be between 0 and 1"))
        new(proportion_sampled, false_positive_rate, false_negative_rate)
    end
end
