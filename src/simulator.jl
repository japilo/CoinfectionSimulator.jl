"""
Simulator module for multi-strain coinfection dynamics.
"""

# Import specific functions we need from StatsBase and Distributions
using StatsBase: sample
using Distributions: Binomial, Poisson

"""
    simulate(initial_population::Population, params::SimulationParameters) -> Vector{Population}

    Simulates multiple pathogen strains spreading through a host population with coinfection dynamics among the strains. Supports both density-dependent and frequency-dependent transmission.

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
    - `transmission_type::Symbol`: Type of transmission (:density or :frequency, defaults to :frequency)

# Returns
Vector of population states at each time step (Vector{Population})
"""
function simulate(initial_population::Population, params::SimulationParameters)::Vector{Population}
    # Set up results collection
    n_steps = params.time_steps
    result_pop = Vector{Population}(undef, n_steps)
    working_pop = copy_population(initial_population)
    result_pop[1] = copy_population(working_pop)
    mortality_list = Set{Int}()
    alive_indices = Vector{Int}()

    # Strain introduction timing
    n_strains = isempty(initial_population) ? length(params.models) : size(initial_population[1], 1)
    intro_step = if params.introduction == :simultaneous
        ones(Int, n_strains)
    elseif params.introduction == :random
        rand(1:n_steps, n_strains)
    else # :none
        zeros(Int, n_strains)
    end

    for t in 1:(n_steps-1)
        # Garbage collection
        empty!(mortality_list)

        # Introduce infections
        if any(intro_step .== t)
            introduce_infections!(working_pop, t, intro_step)
        end

        # Breeding
        breeding!(working_pop, params.fecundity, params.age_maturity)

        # Process disease dynamics
        process_disease_dynamics!(working_pop, params, mortality_list)

        # Remove dead individuals
        if !isempty(mortality_list)
            remove_dead_individuals!(working_pop, mortality_list, alive_indices)
        end

        # Age surviving individuals in-place
        age_population!(working_pop)

        # Store snapshot of current state
        result_pop[t+1] = copy_population(working_pop)
    end

    return result_pop
end
