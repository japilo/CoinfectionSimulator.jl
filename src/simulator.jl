"""
Simulator module for multi-strain coinfection dynamics.
"""

# Import specific functions we need from StatsBase and Distributions
using StatsBase: sample
using Distributions: Binomial, Poisson

"""
    simulate(initial_population::Population, params::SimulationParameters) -> Vector{Population}

    Simulates multiple pathogen strains spreading through a host population with coinfection dynamics among the strains. Density-dependent transmission is implemented.

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
Vector of population states at each time step (Vector{Population})
"""
function simulate(initial_population::Population, params::SimulationParameters)
    # Initialize results
    n_steps = params.time_steps
    result_pop = Vector{Population}(undef, n_steps)
    result_pop[1] = deepcopy(initial_population)

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
    end

    return result_pop
end
