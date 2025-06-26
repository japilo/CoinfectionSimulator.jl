"""
Virtual ecologist sampling module - simulates imperfect detection in ecological surveys.
"""

"""
    virtual_ecologist_sample(;
        virtual_population,
        proportion_sampled::Float64,
        false_positive_rate::Float64,
        false_negative_rate::Float64
    ) -> Matrix{Bool}

Simulates sampling by a virtual ecologist with imperfect detection capabilities.

# Arguments
- `virtual_population`: Population data across time steps (Vector{Vector{Matrix{Bool}}} or similar).
  Each element represents a time step, containing a vector of individual matrices.
  Each individual matrix has rows for strains and columns for disease states [S,E,I,R].
- `proportion_sampled::Float64`: Proportion of the population sampled at each time step (0-1).
- `false_positive_rate::Float64`: Probability of detecting a strain when it's not present (0-1).
- `false_negative_rate::Float64`: Probability of not detecting a strain when it is present (0-1).

# Returns
- `Matrix{Bool}`: Detection matrix with dimensions (n_timesteps, n_strains).
  Each element indicates whether a strain was detected at that time step.

# Examples
```julia
# Create sample population data
pop = [[rand(Bool, 3, 4) for _ in 1:10] for _ in 1:5]  # 5 timesteps, 10 individuals, 3 strains

# Sample with 50% sampling rate and 10% error rates
detections = virtual_ecologist_sample(
    virtual_population=pop,
    proportion_sampled=0.5,
    false_positive_rate=0.1,
    false_negative_rate=0.1
)
```
"""
function virtual_ecologist_sample(;
	virtual_population,
	proportion_sampled::Float64,
	false_positive_rate::Float64,
	false_negative_rate::Float64,
)
	# Input validation
	@assert length(virtual_population) > 0 "The virtual population must not be empty."
	@assert all(t -> all(m -> size(m, 2) == 4, t), virtual_population) "Each matrix in the virtual population must have 4 columns."
	@assert all(t -> length(unique([size(m, 1) for m in t])) <= 1, virtual_population) "All matrices within each timestep must have the same number of rows."
	@assert proportion_sampled >= 0 && proportion_sampled <= 1 "Proportion sampled must be between 0 and 1."
	@assert false_positive_rate >= 0 && false_positive_rate <= 1 "False positive rate must be between 0 and 1."
	@assert false_negative_rate >= 0 && false_negative_rate <= 1 "False negative rate must be between 0 and 1."

	# Get dimensions
	n_timesteps = length(virtual_population)
	n_strains = isempty(virtual_population[1]) ? 0 : size(virtual_population[1][1], 1)
	
	# Initialize detection matrix
	detect_matrix = falses(n_timesteps, n_strains)

	# Process each timestep
	for t in 1:n_timesteps
		timestep_pop = virtual_population[t]
		
		# Skip empty timesteps
		if isempty(timestep_pop)
			continue
		end
		
		n_individuals = length(timestep_pop)
		individuals_sampled = max(1, round(Int, proportion_sampled * n_individuals))
		
		# Randomly sample individuals
		sampled_indices = sample(1:n_individuals, individuals_sampled; replace=false)
		
		# Extract infection status (column 3) for sampled individuals
		perfect_sample = [timestep_pop[i][:, 3] for i in sampled_indices]
		imperfect_sample = falses(individuals_sampled, n_strains)

		# Apply false positive and false negative rates
		for i in 1:individuals_sampled
			for j in 1:n_strains
				if perfect_sample[i][j] == true
					# Individual is truly infected - apply false negative rate
					imperfect_sample[i, j] = rand() > false_negative_rate
				else
					# Individual is not infected - apply false positive rate
					imperfect_sample[i, j] = rand() < false_positive_rate
				end
			end
		end

		# Determine strain detection at timestep level
		# A strain is detected if at least one sampled individual tests positive
		detect_matrix[t, :] = vec(sum(imperfect_sample, dims=1)) .> 0
	end

	return detect_matrix
end
