"""
Simulator module for multi-strain coinfection dynamics.
"""

# Import specific functions we need from StatsBase and Distributions
using StatsBase: sample
using Distributions: Binomial, Poisson

"""
    simulate(initial_population::Population, params::SimulationParameters) -> Tuple{Vector{Population}, Vector{Vector{Int}}}

    Simulates multiple pathogen strains spreading through a host population with coinfection dynamics among the strains.

# Arguments
- `initial_population::Population`: Initial population of individuals
- `params::SimulationParameters`: Parameters for the simulation. The fields are:
    - `models::Vector{<:DiseaseModel}`: Disease models for each strain
    - `interactions::Matrix{Float64}`: Strain interaction matrix
    - `base_mortality::Float64`: Background mortality rate
    - `fecundity::Float64`: Mean number of offspring per mature individual
    - `age_maturity::Int`: Age (in time steps) at which individuals can reproduce
    - `introduction::Symbol`: How strains are introduced (:simultaneous, :random, or :none)
    - `time_steps::Int`: Number of time steps to simulate

# Returns
A tuple containing:
1. Vector of population states at each time step (Vector{Population})
2. Vector of individual ages at each time step (Vector{Vector{Int}})
"""
function simulate(initial_population::Population, params::SimulationParameters)
    # Initialize results
    n_steps = params.time_steps
    result_pop = Vector{Population}(undef, n_steps)
    result_pop[1] = deepcopy(initial_population)
    result_ages = Vector{Vector{Int}}(undef, n_steps)
    result_ages[1] = [ind.age for ind in initial_population]

    # Strain introduction timing
    n_strains = isempty(initial_population) ? length(params.models) : size(initial_population[1], 1)
    intro_step = if params.introduction == :simultaneous
        ones(Int, n_strains)
    elseif params.introduction == :random
        rand(1:n_steps, n_strains)
    else # :none
        zeros(Int, n_strains)
    end

    # Time loop
    for t in 1:(n_steps-1)
        current_pop = deepcopy(result_pop[t])

        # Introduce infections
        if any(intro_step .== t)
            introduce_infections!(current_pop, t, intro_step)
        end

        # Breeding
        breeding!(current_pop, params.fecundity, params.age_maturity)

        # Process disease dynamics
        mortality_list = Set{Int}()
        process_disease_dynamics!(current_pop, params, mortality_list)

        # Apply mortality
        if !isempty(mortality_list)
            alive_mask = [i ∉ mortality_list for i in 1:length(current_pop)]
            current_pop = Population(current_pop.individuals[alive_mask])
        end

        # Age surviving individuals
        for ind in current_pop
            ind.age += 1
        end

        # Store results
        result_pop[t+1] = current_pop
        result_ages[t+1] = [ind.age for ind in current_pop]
    end

    return (result_pop, result_ages)
end

"""
    introduce_infections!(population::Population, timestep::Int, intro_schedule::Vector{Int})

Introduce infections into the population according to the introduction schedule.
"""
function introduce_infections!(population::Population, timestep::Int, intro_schedule::Vector{Int})
    strains_to_introduce = findall(intro_schedule .== timestep)
    isempty(strains_to_introduce) && return

    n_individuals = length(population)
    n_to_infect = min(length(strains_to_introduce), n_individuals)

    infected_indices = sample(1:n_individuals, n_to_infect, replace=false)

    for (i, idx) in enumerate(infected_indices)
        strain = strains_to_introduce[i]
        # Change state from susceptible to infected
        population[idx][strain, 1] = false
        population[idx][strain, 3] = true
    end
end

"""
    breeding!(population::Population, fecundity::Float64, maturity_age::Int)

Handle reproduction of mature individuals in the population.
"""
function breeding!(population::Population, fecundity::Float64, maturity_age::Int)
    fecundity <= 0 && return

    # Count breeding-age individuals
    mature_count = count(ind -> ind.age >= maturity_age, population)
    mature_count == 0 && return

    # Sample number of births from Poisson distribution
    n_births = rand(Poisson(mature_count * fecundity))
    n_births == 0 && return

    # Create new individuals
    n_strains = size(population[1], 1)
    for _ in 1:n_births
        add_individual!(population, Individual(n_strains, 0))
    end
end

"""
    process_disease_dynamics!(population::Population, params::SimulationParameters, mortality_list::Set{Int})

Process all disease dynamics for one time step.
"""
function process_disease_dynamics!(population::Population, params::SimulationParameters, mortality_list::Set{Int})
    isempty(population) && return

    # For each strain, check if any infections are active
    n_strains = length(params.models)
    active_infections = falses(n_strains)

    for ind in population
        for strain in 1:n_strains
            if ind[strain, 2] || ind[strain, 3]  # Exposed or Infected
                active_infections[strain] = true
                break
            end
        end
    end

    # Get alive individuals (not already marked for death)
    alive = [i ∉ mortality_list for i in 1:length(population)]

    # Process each strain
    for strain in 1:n_strains
        # Base mortality for uninfected individuals
        apply_base_mortality!(population, alive, strain, params.base_mortality, mortality_list)

        # Skip if no active infections for this strain
        !active_infections[strain] && continue

        # Update alive mask
        alive = [i ∉ mortality_list for i in 1:length(population)]

        # Process the appropriate disease model
        process_strain!(population, alive, strain, params.models[strain], params.interactions[strain, :],
            params.base_mortality, mortality_list)
    end
end

"""
    apply_base_mortality!(population::Population, alive::Vector{Bool}, strain::Int,
                         base_mortality::Float64, mortality_list::Set{Int})

Apply background mortality to uninfected individuals.
"""
function apply_base_mortality!(population::Population, alive::Vector{Bool}, strain::Int,
    base_mortality::Float64, mortality_list::Set{Int})
    base_mortality <= 0 && return

    # Find uninfected individuals (in S, E, or R state)
    uninfected_indices = findall(i -> alive[i] && !population[i][strain, 3], 1:length(population))
    isempty(uninfected_indices) && return

    # Apply mortality
    n_deaths = rand(Binomial(length(uninfected_indices), base_mortality))
    n_deaths <= 0 && return

    dead_indices = sample(uninfected_indices, n_deaths, replace=false)
    union!(mortality_list, dead_indices)
end


"""
    process_strain!(population::Population, alive::Vector{Bool}, strain::Int,
                   model::DiseaseModel, interactions::Vector{Float64},
                   base_mortality::Float64, mortality_list::Set{Int})

Process disease dynamics for a specific strain based on its disease model.
"""
function process_strain!(population::Population, alive::Vector{Bool}, strain::Int,
    model::SIModel, interactions::Vector{Float64},
    base_mortality::Float64, mortality_list::Set{Int})
    # Process infection spread
    process_infections!(population, alive, strain, model.transmission, interactions)

    # Process mortality from infection
    process_disease_mortality!(population, alive, strain, base_mortality, model.mortality, mortality_list)
end

function process_strain!(population::Population, alive::Vector{Bool}, strain::Int,
    model::SIRModel, interactions::Vector{Float64},
    base_mortality::Float64, mortality_list::Set{Int})
    # Process infection spread
    process_infections!(population, alive, strain, model.transmission, interactions)

    # Process recovery
    process_recovery!(population, alive, strain, model.recovery)

    # Process mortality from infection
    process_disease_mortality!(population, alive, strain, base_mortality, model.mortality, mortality_list)
end

function process_strain!(population::Population, alive::Vector{Bool}, strain::Int,
    model::SEIRModel, interactions::Vector{Float64},
    base_mortality::Float64, mortality_list::Set{Int})
    # Process exposure
    process_exposures!(population, alive, strain, model.transmission, interactions)

    # Process transition from exposed to infected
    process_latent_infections!(population, alive, strain, model.latency)

    # Process recovery
    process_recovery!(population, alive, strain, model.recovery)

    # Process mortality from infection
    process_disease_mortality!(population, alive, strain, base_mortality, model.mortality, mortality_list)
end

function process_strain!(population::Population, alive::Vector{Bool}, strain::Int,
    model::SEIRSModel, interactions::Vector{Float64},
    base_mortality::Float64, mortality_list::Set{Int})
    # Process exposure
    process_exposures!(population, alive, strain, model.transmission, interactions)

    # Process transition from exposed to infected
    process_latent_infections!(population, alive, strain, model.latency)

    # Process recovery
    process_recovery!(population, alive, strain, model.recovery)

    # Process immunity loss
    process_immunity_loss!(population, alive, strain, model.immunity_loss)

    # Process mortality from infection
    process_disease_mortality!(population, alive, strain, base_mortality, model.mortality, mortality_list)
end

"""
    process_infections!(population::Population, alive::Vector{Bool}, strain::Int,
                       transmission_rate::Float64, interactions::Vector{Float64})

Process direct infections for SI/SIR models (susceptible to infected).
"""
function process_infections!(population::Population, alive::Vector{Bool}, strain::Int,
    transmission_rate::Float64, interactions::Vector{Float64})
    transmission_rate <= 0 && return

    # Find infected and susceptible individuals
    infected_indices = findall(i -> alive[i] && population[i][strain, 3], 1:length(population))
    susceptible_indices = findall(i -> alive[i] && population[i][strain, 1], 1:length(population))

    isempty(infected_indices) && return
    isempty(susceptible_indices) && return

    # Calculate infection probability for each susceptible individual
    n_infected = length(infected_indices)

    for s_idx in susceptible_indices
        # Get the baseline infection pressure
        infection_pressure = transmission_rate * n_infected / length(population)

        # Adjust for coinfection interactions
        for other_strain in 1:length(interactions)
            other_strain == strain && continue

            if population[s_idx][other_strain, 3]  # Individual is infected with other strain
                infection_pressure *= interactions[other_strain]
            end
        end

        # Apply infection
        if rand() < min(1.0, infection_pressure)
            population[s_idx][strain, 1] = false  # No longer susceptible
            population[s_idx][strain, 3] = true   # Now infected
        end
    end
end

"""
    process_exposures!(population::Population, alive::Vector{Bool}, strain::Int,
                      transmission_rate::Float64, interactions::Vector{Float64})

Process exposures for SEIR/SEIRS models (susceptible to exposed).
"""
function process_exposures!(population::Population, alive::Vector{Bool}, strain::Int,
    transmission_rate::Float64, interactions::Vector{Float64})
    transmission_rate <= 0 && return

    # Find infected and susceptible individuals
    infected_indices = findall(i -> alive[i] && population[i][strain, 3], 1:length(population))
    susceptible_indices = findall(i -> alive[i] && population[i][strain, 1], 1:length(population))

    isempty(infected_indices) && return
    isempty(susceptible_indices) && return

    # Calculate infection probability for each susceptible individual
    n_infected = length(infected_indices)

    for s_idx in susceptible_indices
        # Get the baseline infection pressure
        infection_pressure = transmission_rate * n_infected / length(population)

        # Adjust for coinfection interactions
        for other_strain in 1:length(interactions)
            other_strain == strain && continue

            if population[s_idx][other_strain, 3]  # Individual is infected with other strain
                infection_pressure *= interactions[other_strain]
            end
        end

        # Apply exposure
        if rand() < min(1.0, infection_pressure)
            population[s_idx][strain, 1] = false  # No longer susceptible
            population[s_idx][strain, 2] = true   # Now exposed
        end
    end
end

"""
    process_latent_infections!(population::Population, alive::Vector{Bool}, strain::Int, latency::Int)

Process transitions from exposed to infected based on latency period.
"""
function process_latent_infections!(population::Population, alive::Vector{Bool}, strain::Int, latency::Int)
    # This is a simple implementation that transitions exposed to infected with probability 1/latency
    # A more sophisticated version might track exposure time for each individual

    exposed_indices = findall(i -> alive[i] && population[i][strain, 2], 1:length(population))
    isempty(exposed_indices) && return

    transition_prob = 1.0 / latency

    for e_idx in exposed_indices
        if rand() < transition_prob
            population[e_idx][strain, 2] = false  # No longer exposed
            population[e_idx][strain, 3] = true   # Now infected
        end
    end
end

"""
    process_recovery!(population::Population, alive::Vector{Bool}, strain::Int, recovery_rate::Float64)

Process recovery from infection.
"""
function process_recovery!(population::Population, alive::Vector{Bool}, strain::Int, recovery_rate::Float64)
    recovery_rate <= 0 && return

    infected_indices = findall(i -> alive[i] && population[i][strain, 3], 1:length(population))
    isempty(infected_indices) && return

    for i_idx in infected_indices
        if rand() < recovery_rate
            population[i_idx][strain, 3] = false  # No longer infected
            population[i_idx][strain, 4] = true   # Now recovered
        end
    end
end

"""
    process_immunity_loss!(population::Population, alive::Vector{Bool}, strain::Int, immunity_loss_rate::Float64)

Process loss of immunity after recovery.
"""
function process_immunity_loss!(population::Population, alive::Vector{Bool}, strain::Int, immunity_loss_rate::Float64)
    immunity_loss_rate <= 0 && return

    recovered_indices = findall(i -> alive[i] && population[i][strain, 4], 1:length(population))
    isempty(recovered_indices) && return

    for r_idx in recovered_indices
        if rand() < immunity_loss_rate
            population[r_idx][strain, 4] = false  # No longer recovered
            population[r_idx][strain, 1] = true   # Now susceptible again
        end
    end
end

"""
    process_disease_mortality!(population::Population, alive::Vector{Bool}, strain::Int,
                              base_mortality::Float64, disease_mortality::Float64, mortality_list::Set{Int})

Process additional mortality due to infection.
"""
function process_disease_mortality!(population::Population, alive::Vector{Bool}, strain::Int,
    base_mortality::Float64, disease_mortality::Float64, mortality_list::Set{Int})
    disease_mortality <= 0 && return

    infected_indices = findall(i -> alive[i] && population[i][strain, 3], 1:length(population))
    isempty(infected_indices) && return

    # Calculate total mortality probability (base + disease)
    total_mortality = base_mortality + disease_mortality

    # Apply mortality
    n_deaths = rand(Binomial(length(infected_indices), total_mortality))
    n_deaths <= 0 && return

    dead_indices = sample(infected_indices, n_deaths, replace=false)
    union!(mortality_list, dead_indices)
end

# Backwards compatibility function
"""
    coinfection_simulator(;
        initial_pop::Vector{<:AbstractMatrix{Bool}},
        ages::Vector{Int},
        interactions::Matrix{Float64},
        disease_type::Vector{String},
        base_mortality::Float64,
        disease_mortality::Vector{Float64},
        fecundity::Float64,
        transmission::Vector{Float64},
        time_steps::Int,
        age_maturity::Int,
        introduction::String = "simultaneous",
        latency::Union{Vector{Int}, Nothing} = nothing,
        recovery::Union{Vector{Float64}, Nothing} = nothing,
        immunity_loss::Union{Vector{Float64}, Nothing} = nothing
    ) -> Tuple{Vector{Vector{<:AbstractMatrix{Bool}}}, Vector{Vector{Int}}}

Legacy interface for the coinfection simulator. Converts parameters to the new type system
and calls the `simulate` function.
"""
function coinfection_simulator(;
    initial_pop::Vector{<:AbstractMatrix{Bool}},
    ages::Vector{Int},
    interactions::Matrix{Float64},
    disease_type::Vector{String},
    base_mortality::Float64,
    disease_mortality::Vector{Float64},
    fecundity::Float64,
    transmission::Vector{Float64},
    time_steps::Int,
    age_maturity::Int,
    introduction::String="simultaneous",
    latency::Union{Vector{Int},Nothing}=nothing,
    recovery::Union{Vector{Float64},Nothing}=nothing,
    immunity_loss::Union{Vector{Float64},Nothing}=nothing
)
    # Convert string introduction type to symbol
    intro_sym = if introduction == "simultaneous"
        :simultaneous
    elseif introduction == "random"
        :random
    elseif introduction == "none"
        :none
    else
        throw(ArgumentError("Invalid introduction type: $introduction"))
    end

    # Convert disease models
    models = Vector{DiseaseModel}(undef, length(disease_type))
    for i in 1:length(disease_type)
        if disease_type[i] == "si"
            models[i] = SIModel(transmission[i], disease_mortality[i])
        elseif disease_type[i] == "sir"
            recovery === nothing && throw(ArgumentError("Recovery rates must be provided for SIR models"))
            models[i] = SIRModel(transmission[i], disease_mortality[i], recovery[i])
        elseif disease_type[i] == "seir"
            recovery === nothing && throw(ArgumentError("Recovery rates must be provided for SEIR models"))
            latency === nothing && throw(ArgumentError("Latency periods must be provided for SEIR models"))
            models[i] = SEIRModel(transmission[i], disease_mortality[i], recovery[i], latency[i])
        elseif disease_type[i] == "seirs"
            recovery === nothing && throw(ArgumentError("Recovery rates must be provided for SEIRS models"))
            latency === nothing && throw(ArgumentError("Latency periods must be provided for SEIRS models"))
            immunity_loss === nothing && throw(ArgumentError("Immunity loss rates must be provided for SEIRS models"))
            models[i] = SEIRSModel(transmission[i], disease_mortality[i], recovery[i], latency[i], immunity_loss[i])
        else
            throw(ArgumentError("Invalid disease type: $(disease_type[i])"))
        end
    end

    # Create simulation parameters
    params = SimulationParameters(
        models,
        interactions,
        base_mortality,
        fecundity,
        age_maturity,
        intro_sym,
        time_steps
    )

    # Create initial population
    initial_population = Population(initial_pop, ages)

    # Run simulation
    result_pop, result_ages = simulate(initial_population, params)

    # Convert results back to legacy format
    legacy_pop = Vector{Vector{eltype(initial_pop)}}(undef, time_steps)
    for t in 1:time_steps
        legacy_pop[t] = [ind.state for ind in result_pop[t]]
    end

    return (legacy_pop, result_ages)
end
