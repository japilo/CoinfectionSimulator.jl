"""
Helper functions for the simulator.
"""

# Import specific functions we need from StatsBase and Distributions
using StatsBase: sample
using Distributions: Binomial, Poisson

"""
    introduce_infections!(population::Population, timestep::Int, intro_schedule::Vector{Int})

Introduce infections into the population according to the introduction schedule. Individuals are randomly chosen and their state is changed from susceptible to infected.
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

Handle reproduction of mature individuals in the population. All individuals at or greater than the age of maturity are able to produce offspring. New individuals are generated via a Poisson distribution and added to the population as susceptibles.
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

Process all disease dynamics for one time step. This function applies base mortality to all individuals, checks which strains are active in the population, and handles disease dynamics for those strains.
"""
function process_disease_dynamics!(population::Population, params::SimulationParameters, mortality_list::Set{Int})
    isempty(population) && return

    # For each strain, check if any infections are active or if it's a SEIRS model with recovered individuals
    n_strains = length(params.models)
    active_infections = falses(n_strains)
    has_recovered = falses(n_strains)

    for ind in population
        for strain in 1:n_strains
            if ind[strain, 2] || ind[strain, 3]  # Exposed or Infected
                active_infections[strain] = true
            elseif ind[strain, 4] && isa(params.models[strain], SEIRSModel)  # Recovered and SEIRS model
                has_recovered[strain] = true
            end
        end
    end

    # Apply base mortality
    alive = [i ∉ mortality_list for i in 1:length(population)]
    apply_base_mortality!(population, alive, params.base_mortality, mortality_list)

    # Update alive status after base mortality
    alive = [i ∉ mortality_list for i in 1:length(population)]

    # Process each strain
    for strain in 1:n_strains
        # Skip if no active infections AND not a SEIRS model with recovered individuals
        if !active_infections[strain] && (!has_recovered[strain] || !isa(params.models[strain], SEIRSModel))
            continue
        end

        # Update alive mask
        alive = [i ∉ mortality_list for i in 1:length(population)]

        # Process the appropriate disease model
        process_strain!(population, alive, strain, params.models[strain], params.interactions[strain, :],
            mortality_list)
    end
end

"""
    apply_base_mortality!(population::Population, alive::Vector{Bool}, base_mortality::Float64, mortality_list::Set{Int})

Apply background mortality to all living individuals once per time step.
"""
function apply_base_mortality!(population::Population, alive::Vector{Bool},
    base_mortality::Float64, mortality_list::Set{Int})
    base_mortality <= 0 && return

    # Find all living individuals
    alive_indices = findall(i -> alive[i], 1:length(population))
    isempty(alive_indices) && return

    # Apply mortality
    n_deaths = rand(Binomial(length(alive_indices), base_mortality))
    n_deaths <= 0 && return

    dead_indices = sample(alive_indices, n_deaths, replace=false)
    union!(mortality_list, dead_indices)
end

"""
    process_strain!(population::Population, alive::Vector{Bool}, strain::Int,
                   model::DiseaseModel, interactions::Vector{Float64},
                   base_mortality::Float64, mortality_list::Set{Int})

Process disease dynamics for a specific strain based on its disease model. This function uses multiple dispatch to apply the appropriate disease processes to each model.
"""
function process_strain!(population::Population, alive::Vector{Bool}, strain::Int,
    model::SIModel, interactions::Vector{Float64}, mortality_list::Set{Int})
    # Process infection spread
    process_infections!(population, alive, strain, model.transmission, interactions)

    # Process mortality from infection (base mortality handled separately)
    process_disease_mortality!(population, alive, strain, model.mortality, mortality_list)
end

function process_strain!(population::Population, alive::Vector{Bool}, strain::Int,
    model::SIRModel, interactions::Vector{Float64}, mortality_list::Set{Int})
    # Process infection spread
    process_infections!(population, alive, strain, model.transmission, interactions)

    # Process recovery
    process_recovery!(population, alive, strain, model.recovery)

    # Process mortality from infection (base mortality handled separately)
    process_disease_mortality!(population, alive, strain, model.mortality, mortality_list)
end

function process_strain!(population::Population, alive::Vector{Bool}, strain::Int,
    model::SEIRModel, interactions::Vector{Float64}, mortality_list::Set{Int})
    # Process exposure
    process_exposures!(population, alive, strain, model.transmission, interactions)

    # Process transition from exposed to infected
    process_latent_infections!(population, alive, strain, model.latency)

    # Process recovery
    process_recovery!(population, alive, strain, model.recovery)

    # Process mortality from infection (base mortality handled separately)
    process_disease_mortality!(population, alive, strain, model.mortality, mortality_list)
end

function process_strain!(population::Population, alive::Vector{Bool}, strain::Int,
    model::SEIRSModel, interactions::Vector{Float64}, mortality_list::Set{Int})
    # Process exposure
    process_exposures!(population, alive, strain, model.transmission, interactions)

    # Process transition from exposed to infected
    process_latent_infections!(population, alive, strain, model.latency)

    # Process recovery
    process_recovery!(population, alive, strain, model.recovery)

    # Process immunity loss
    process_immunity_loss!(population, alive, strain, model.immunity_loss)

    # Process mortality from infection (base mortality handled separately)
    process_disease_mortality!(population, alive, strain, model.mortality, mortality_list)
end

"""
    process_infections!(population::Population, alive::Vector{Bool}, strain::Int,
                       transmission_rate::Float64, interactions::Vector{Float64})

Process direct infections for SI/SIR models (susceptible to infected, with no latent period). Interactions with other strains in potential new hosts affect the transmission rate.
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

Process exposures for SEIR/SEIRS models (susceptible to exposed). If a potential new host is infected with another strain, interactions between strains affect the infection pressure.
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

Process transitions from exposed to infected with a probability of 1/latency.
"""
function process_latent_infections!(population::Population, alive::Vector{Bool}, strain::Int, latency::Int)
    # This is a simple implementation that transitions exposed to infected with probability 1/latency

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
                              disease_mortality::Float64, mortality_list::Set{Int})

Process additional mortality due to infection.
"""
function process_disease_mortality!(population::Population, alive::Vector{Bool}, strain::Int,
    disease_mortality::Float64, mortality_list::Set{Int})
    disease_mortality <= 0 && return

    infected_indices = findall(i -> alive[i] && population[i][strain, 3], 1:length(population))
    isempty(infected_indices) && return

    # Apply disease mortality
    n_deaths = rand(Binomial(length(infected_indices), disease_mortality))
    n_deaths <= 0 && return

    dead_indices = sample(infected_indices, n_deaths, replace=false)
    union!(mortality_list, dead_indices)
end
