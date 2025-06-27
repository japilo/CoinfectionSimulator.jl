"""
Utility functions for disease state handling in the coinfection simulator.
"""

using Random
using Distributions
using StatsBase

"""
    infect(susceptibles, infecteds, interactions, beta, strain)

Core infection function that determines which susceptible individuals become infected
based on transmission probability and strain interactions.
"""
function infect(
	susceptibles::Vector{BitMatrix},
	infecteds::Vector{BitMatrix},
	interactions::Vector{Float64},
	beta::Float64,
	strain::Int,
)
	# Convert susceptibles to individual/strain matrix
	strain_matrix = map(m -> sum.(eachrow(m[:, [2, 3]])), susceptibles)
	strain_matrix = vcat(strain_matrix'...)

	indiv = [strain_matrix[i, :] for i in axes(strain_matrix)[1]]
	total_strains = sum.(indiv)

	# Infect individuals
	indiv_new = Vector{Int}(undef, length(indiv))
	for i in eachindex(indiv)
		indiv_current = indiv[i]
		int_indiv = interactions .* indiv_current

		if total_strains[i] > 1
			prob = prod(int_indiv) * beta
			indiv_current[strain] = max(indiv_current[strain], rand(Binomial(length(infecteds), prob)))
		elseif total_strains[i] == 1
			prob = sum(int_indiv) * beta
			indiv_current[strain] = max(indiv_current[strain], rand(Binomial(length(infecteds), prob)))
		else
			indiv_current[strain] = rand(Binomial(length(infecteds), beta))
		end

		indiv_new[i] = indiv_current[strain]
	end

	return indiv_new .> 0
end

"""
    handle_si_disease(current_pop, alive, strain, transmission, interactions, base_mortality, disease_mortality, dead_indices)

Handle SI (Susceptible-Infected) disease dynamics for a single strain.
"""
function handle_si_disease(
	current_pop::Vector{BitMatrix},
	alive::Vector{Bool},
	strain::Int,
	transmission::Float64,
	interactions::Vector{Float64},
	base_mortality::Float64,
	disease_mortality::Float64,
	dead_indices::Set{Int},
)
	# Infection
	handle_infection(current_pop, alive, strain, transmission, interactions)

	# Death of infected
	handle_infected_death(current_pop, strain, base_mortality, disease_mortality, dead_indices)
end

"""
    handle_sir_disease(current_pop, alive, strain, transmission, interactions, base_mortality, disease_mortality, recovery, dead_indices)

Handle SIR (Susceptible-Infected-Recovered) disease dynamics for a single strain.
"""
function handle_sir_disease(
	current_pop::Vector{BitMatrix},
	alive::Vector{Bool},
	strain::Int,
	transmission::Float64,
	interactions::Vector{Float64},
	base_mortality::Float64,
	disease_mortality::Float64,
	recovery::Float64,
	dead_indices::Set{Int},
)
	# Infection
	handle_infection(current_pop, alive, strain, transmission, interactions)

	# Death of infected
	handle_infected_death(current_pop, strain, base_mortality, disease_mortality, dead_indices)

	# Recovery
	handle_recovery(current_pop, alive, strain, recovery)
end

"""
    handle_seir_disease(current_pop, alive, strain, transmission, interactions, base_mortality, disease_mortality, recovery, latency, dead_indices)

Handle SEIR (Susceptible-Exposed-Infected-Recovered) disease dynamics for a single strain.
"""
function handle_seir_disease(
	current_pop::Vector{BitMatrix},
	alive::Vector{Bool},
	strain::Int,
	transmission::Float64,
	interactions::Vector{Float64},
	base_mortality::Float64,
	disease_mortality::Float64,
	recovery::Float64,
	latency::Int,
	dead_indices::Set{Int},
)
	# Exposure
	handle_exposure(current_pop, alive, strain, transmission, interactions)

	# Infection from exposed
	handle_exposed_infection(current_pop, alive, strain, latency)

	# Death of infected
	handle_infected_death(current_pop, strain, base_mortality, disease_mortality, dead_indices)

	# Recovery
	handle_recovery(current_pop, alive, strain, recovery)
end

"""
    handle_seirs_disease(current_pop, alive, strain, transmission, interactions, base_mortality, disease_mortality, recovery, latency, immunity_loss, dead_indices)

Handle SEIRS (Susceptible-Exposed-Infected-Recovered-Susceptible) disease dynamics for a single strain.
"""
function handle_seirs_disease(
	current_pop::Vector{BitMatrix},
	alive::Vector{Bool},
	strain::Int,
	transmission::Float64,
	interactions::Vector{Float64},
	base_mortality::Float64,
	disease_mortality::Float64,
	recovery::Float64,
	latency::Int,
	immunity_loss::Float64,
	dead_indices::Set{Int},
)
	# SEIR processes
	handle_seir_disease(current_pop, alive, strain, transmission, interactions,
		base_mortality, disease_mortality, recovery, latency, dead_indices)

	# Loss of immunity
	handle_immunity_loss(current_pop, alive, strain, immunity_loss)
end

"""
    handle_infection(current_pop, alive, strain, transmission, interactions)

Handle the infection process for susceptible individuals.
"""
function handle_infection(
	current_pop::Vector{BitMatrix},
	alive::Vector{Bool},
	strain::Int,
	transmission::Float64,
	interactions::Vector{Float64},
)
	susceptible_indices = findall(i -> alive[i] && current_pop[i][strain, 1], 1:length(current_pop))
	infected_indices = findall(i -> alive[i] && current_pop[i][strain, 3], 1:length(current_pop))

	if !isempty(susceptible_indices) && !isempty(infected_indices)
		susceptible_pop = current_pop[susceptible_indices]
		infected_pop = current_pop[infected_indices]
		
		new_infections = infect(susceptible_pop, infected_pop, interactions, transmission, strain)
		
		for (i, is_infected) in enumerate(new_infections)
			if is_infected
				idx = susceptible_indices[i]
				current_pop[idx][strain, 1] = false
				current_pop[idx][strain, 3] = true
			end
		end
	end
end

"""
    handle_exposure(current_pop, alive, strain, transmission, interactions)

Handle the exposure process for susceptible individuals in SEIR models.
"""
function handle_exposure(
	current_pop::Vector{BitMatrix},
	alive::Vector{Bool},
	strain::Int,
	transmission::Float64,
	interactions::Vector{Float64},
)
	susceptible_indices = findall(i -> alive[i] && current_pop[i][strain, 1], 1:length(current_pop))
	infected_indices = findall(i -> alive[i] && current_pop[i][strain, 3], 1:length(current_pop))

	if !isempty(susceptible_indices) && !isempty(infected_indices)
		susceptible_pop = current_pop[susceptible_indices]
		infected_pop = current_pop[infected_indices]
		
		new_exposures = infect(susceptible_pop, infected_pop, interactions, transmission, strain)
		
		for (i, is_exposed) in enumerate(new_exposures)
			if is_exposed
				idx = susceptible_indices[i]
				current_pop[idx][strain, 1] = false
				current_pop[idx][strain, 2] = true
			end
		end
	end
end

"""
    handle_exposed_infection(current_pop, alive, strain, latency)

Handle the transition from exposed to infected state.
"""
function handle_exposed_infection(
	current_pop::Vector{BitMatrix},
	alive::Vector{Bool},
	strain::Int,
	latency::Int,
)
	exposed_indices = findall(i -> alive[i] && current_pop[i][strain, 2], 1:length(current_pop))
	
	if !isempty(exposed_indices)
		infection_prob = 1.0 / latency
		n_infections = rand(Binomial(length(exposed_indices), infection_prob))
		
		if n_infections > 0
			infected_indices = sample(exposed_indices, n_infections; replace = false)
			for idx in infected_indices
				current_pop[idx][strain, 2] = false
				current_pop[idx][strain, 3] = true
			end
		end
	end
end

"""
    handle_recovery(current_pop, alive, strain, recovery)

Handle recovery from infected to recovered state.
"""
function handle_recovery(
	current_pop::Vector{BitMatrix},
	alive::Vector{Bool},
	strain::Int,
	recovery::Float64,
)
	infected_indices = findall(i -> alive[i] && current_pop[i][strain, 3], 1:length(current_pop))
	
	if !isempty(infected_indices)
		n_recoveries = rand(Binomial(length(infected_indices), recovery))
		
		if n_recoveries > 0
			recovered_indices = sample(infected_indices, n_recoveries; replace = false)
			for idx in recovered_indices
				current_pop[idx][strain, 3] = false
				current_pop[idx][strain, 4] = true
			end
		end
	end
end

"""
    handle_immunity_loss(current_pop, alive, strain, immunity_loss)

Handle loss of immunity (transition from recovered to susceptible).
"""
function handle_immunity_loss(
	current_pop::Vector{BitMatrix},
	alive::Vector{Bool},
	strain::Int,
	immunity_loss::Float64,
)
	recovered_indices = findall(i -> alive[i] && current_pop[i][strain, 4], 1:length(current_pop))
	
	if !isempty(recovered_indices)
		n_immunity_loss = rand(Binomial(length(recovered_indices), immunity_loss))
		
		if n_immunity_loss > 0
			susceptible_indices = sample(recovered_indices, n_immunity_loss; replace = false)
			for idx in susceptible_indices
				current_pop[idx][strain, 4] = false
				current_pop[idx][strain, 1] = true
			end
		end
	end
end

"""
    handle_infected_death(current_pop, strain, base_mortality, disease_mortality, dead_indices)

Handle mortality of infected individuals.
"""
function handle_infected_death(
	current_pop::Vector{BitMatrix},
	strain::Int,
	base_mortality::Float64,
	disease_mortality::Float64,
	dead_indices::Set{Int},
)
	infected_indices = findall(i -> i ∉ dead_indices && current_pop[i][strain, 3], 1:length(current_pop))
	
	if !isempty(infected_indices)
		total_mortality = base_mortality + disease_mortality
		n_deaths = rand(Binomial(length(infected_indices), total_mortality))
		
		if n_deaths > 0
			dead_individuals = sample(infected_indices, n_deaths; replace = false)
			union!(dead_indices, dead_individuals)
		end
	end
end

# Vector{Matrix{Bool}} overloads for simulator compatibility

function handle_si_disease(
	current_pop::Vector{Matrix{Bool}},
	alive::Vector{Bool},
	strain::Int,
	transmission::Float64,
	interactions::Vector{Float64},
	base_mortality::Float64,
	disease_mortality::Float64,
	dead_indices::Set{Int},
)
	# Infection
	handle_infection(current_pop, alive, strain, transmission, interactions)

	# Death of infected
	handle_infected_death(current_pop, strain, base_mortality, disease_mortality, dead_indices)
end

function handle_sir_disease(
	current_pop::Vector{Matrix{Bool}},
	alive::Vector{Bool},
	strain::Int,
	transmission::Float64,
	interactions::Vector{Float64},
	base_mortality::Float64,
	disease_mortality::Float64,
	recovery::Float64,
	dead_indices::Set{Int},
)
	# Infection
	handle_infection(current_pop, alive, strain, transmission, interactions)

	# Death of infected
	handle_infected_death(current_pop, strain, base_mortality, disease_mortality, dead_indices)

	# Recovery
	handle_recovery(current_pop, alive, strain, recovery)
end

function handle_seir_disease(
	current_pop::Vector{Matrix{Bool}},
	alive::Vector{Bool},
	strain::Int,
	transmission::Float64,
	interactions::Vector{Float64},
	base_mortality::Float64,
	disease_mortality::Float64,
	recovery::Float64,
	latency::Int,
	dead_indices::Set{Int},
)
	# Exposure
	handle_exposure(current_pop, alive, strain, transmission, interactions)

	# Infection from exposed
	handle_exposed_infection(current_pop, alive, strain, latency)

	# Death of infected
	handle_infected_death(current_pop, strain, base_mortality, disease_mortality, dead_indices)

	# Recovery
	handle_recovery(current_pop, alive, strain, recovery)
end

function handle_seirs_disease(
	current_pop::Vector{Matrix{Bool}},
	alive::Vector{Bool},
	strain::Int,
	transmission::Float64,
	interactions::Vector{Float64},
	base_mortality::Float64,
	disease_mortality::Float64,
	recovery::Float64,
	latency::Int,
	immunity_loss::Float64,
	dead_indices::Set{Int},
)
	# SEIR processes
	handle_seir_disease(current_pop, alive, strain, transmission, interactions,
		base_mortality, disease_mortality, recovery, latency, dead_indices)

	# Loss of immunity
	handle_immunity_loss(current_pop, alive, strain, immunity_loss)
end

function handle_infection(
	current_pop::Vector{Matrix{Bool}},
	alive::Vector{Bool},
	strain::Int,
	transmission::Float64,
	interactions::Vector{Float64},
)
	susceptible_indices = findall(i -> alive[i] && current_pop[i][strain, 1], 1:length(current_pop))
	infected_indices = findall(i -> alive[i] && current_pop[i][strain, 3], 1:length(current_pop))

	if !isempty(susceptible_indices) && !isempty(infected_indices)
		susceptible_pop = current_pop[susceptible_indices]
		infected_pop = current_pop[infected_indices]
		
		new_infections = infect(susceptible_pop, infected_pop, interactions, transmission, strain)
		
		for (i, is_infected) in enumerate(new_infections)
			if is_infected
				idx = susceptible_indices[i]
				current_pop[idx][strain, 1] = false
				current_pop[idx][strain, 3] = true
			end
		end
	end
end

function handle_exposure(
	current_pop::Vector{Matrix{Bool}},
	alive::Vector{Bool},
	strain::Int,
	transmission::Float64,
	interactions::Vector{Float64},
)
	susceptible_indices = findall(i -> alive[i] && current_pop[i][strain, 1], 1:length(current_pop))
	infected_indices = findall(i -> alive[i] && current_pop[i][strain, 3], 1:length(current_pop))

	if !isempty(susceptible_indices) && !isempty(infected_indices)
		susceptible_pop = current_pop[susceptible_indices]
		infected_pop = current_pop[infected_indices]
		
		new_exposures = infect(susceptible_pop, infected_pop, interactions, transmission, strain)
		
		for (i, is_exposed) in enumerate(new_exposures)
			if is_exposed
				idx = susceptible_indices[i]
				current_pop[idx][strain, 1] = false
				current_pop[idx][strain, 2] = true
			end
		end
	end
end

function handle_exposed_infection(
	current_pop::Vector{Matrix{Bool}},
	alive::Vector{Bool},
	strain::Int,
	latency::Int,
)
	exposed_indices = findall(i -> alive[i] && current_pop[i][strain, 2], 1:length(current_pop))
	
	if !isempty(exposed_indices)
		infection_prob = 1.0 / latency
		n_infections = rand(Binomial(length(exposed_indices), infection_prob))
		
		if n_infections > 0
			infected_indices = sample(exposed_indices, n_infections; replace = false)
			for idx in infected_indices
				current_pop[idx][strain, 2] = false
				current_pop[idx][strain, 3] = true
			end
		end
	end
end

function handle_recovery(
	current_pop::Vector{Matrix{Bool}},
	alive::Vector{Bool},
	strain::Int,
	recovery::Float64,
)
	infected_indices = findall(i -> alive[i] && current_pop[i][strain, 3], 1:length(current_pop))
	
	if !isempty(infected_indices)
		n_recoveries = rand(Binomial(length(infected_indices), recovery))
		
		if n_recoveries > 0
			recovered_indices = sample(infected_indices, n_recoveries; replace = false)
			for idx in recovered_indices
				current_pop[idx][strain, 3] = false
				current_pop[idx][strain, 4] = true
			end
		end
	end
end

function handle_immunity_loss(
	current_pop::Vector{Matrix{Bool}},
	alive::Vector{Bool},
	strain::Int,
	immunity_loss::Float64,
)
	recovered_indices = findall(i -> alive[i] && current_pop[i][strain, 4], 1:length(current_pop))
	
	if !isempty(recovered_indices)
		n_immunity_loss = rand(Binomial(length(recovered_indices), immunity_loss))
		
		if n_immunity_loss > 0
			susceptible_indices = sample(recovered_indices, n_immunity_loss; replace = false)
			for idx in susceptible_indices
				current_pop[idx][strain, 4] = false
				current_pop[idx][strain, 1] = true
			end
		end
	end
end

function handle_infected_death(
	current_pop::Vector{Matrix{Bool}},
	strain::Int,
	base_mortality::Float64,
	disease_mortality::Float64,
	dead_indices::Set{Int},
)
	infected_indices = findall(i -> i ∉ dead_indices && current_pop[i][strain, 3], 1:length(current_pop))
	
	if !isempty(infected_indices)
		total_mortality = base_mortality + disease_mortality
		n_deaths = rand(Binomial(length(infected_indices), total_mortality))
		
		if n_deaths > 0
			dead_individuals = sample(infected_indices, n_deaths; replace = false)
			union!(dead_indices, dead_individuals)
		end
	end
end

function infect(
	susceptibles::Vector{Matrix{Bool}},
	infecteds::Vector{Matrix{Bool}},
	interactions::Vector{Float64},
	beta::Float64,
	strain::Int,
)
	# Convert susceptibles to individual/strain matrix
	strain_matrix = map(m -> sum.(eachrow(m[:, [2, 3]])), susceptibles)
	strain_matrix = vcat(strain_matrix'...)

	indiv = [strain_matrix[i, :] for i in axes(strain_matrix)[1]]
	total_strains = sum.(indiv)

	# Infect individuals
	indiv_new = Vector{Int}(undef, length(indiv))
	for i in eachindex(indiv)
		indiv_current = indiv[i]
		int_indiv = interactions .* indiv_current

		if total_strains[i] > 1
			prob = prod(int_indiv) * beta
			indiv_current[strain] = max(indiv_current[strain], rand(Binomial(length(infecteds), prob)))
		elseif total_strains[i] == 1
			prob = sum(int_indiv) * beta
			indiv_current[strain] = max(indiv_current[strain], rand(Binomial(length(infecteds), prob)))
		else
			indiv_current[strain] = rand(Binomial(length(infecteds), beta))
		end

		indiv_new[i] = indiv_current[strain]
	end

	return indiv_new .> 0
end
