"""
Virtual ecologist sampling module - simulates imperfect detection in ecological surveys.
"""

"""
    sample_populations(populations::Vector{Population}, params::SamplingParameters) -> Matrix{Bool}

Sample a series of populations using a virtual ecologist with imperfect detection capabilities.

# Arguments
- `populations::Vector{Population}`: Population data across time steps.
- `params::SamplingParameters`: Parameters controlling the sampling process. The fields are:
    - `proportion_sampled::Float64`: Proportion of the population sampled
    - `false_positive_rate::Float64`: Probability of false positive detection
    - `false_negative_rate::Float64`: Probability of false negative detection

# Returns
- `Matrix{Bool}`: Detection matrix with dimensions (n_timesteps, n_strains).
  Each element indicates whether a strain was detected at that time step.

# Examples
```julia
# Create sample population where all individuals are infected with all strains
pop = Population([Individual(BitMatrix([false false true false]), 20) for _ in 1:100])

# Sample with 50% sampling rate and 10% error rates
detections = sample_populations(
    pop,
    SamplingParameters(0.5, 0.1, 0.1)
)
```
"""
function sample_populations(populations::Vector{Population}, params::SamplingParameters)
    isempty(populations) && return Matrix{Bool}(undef, 0, 0)

    # Get dimensions
    n_timesteps = length(populations)
    n_strains = isempty(populations[1]) ? 0 : size(populations[1][1], 1)

    # Initialize detection matrix
    detect_matrix = falses(n_timesteps, n_strains)

    # Process each timestep
    for t in 1:n_timesteps
        current_pop = populations[t]

        # Skip empty timesteps
        isempty(current_pop) && continue

        n_individuals = length(current_pop)
        individuals_sampled = max(1, round(Int, params.proportion_sampled * n_individuals))
        individuals_sampled = min(individuals_sampled, n_individuals)

        # Randomly sample individuals
        sampled_indices = StatsBase.sample(1:n_individuals, individuals_sampled, replace=false)

        # Apply detection process to sampled individuals
        strain_detected = falses(n_strains)

        for idx in sampled_indices
            for strain in 1:n_strains
                # Check if the individual is infected with this strain
                is_infected = current_pop[idx][strain, 3]

                detected = if is_infected
                    # True positive with probability (1 - false_negative_rate)
                    rand() > params.false_negative_rate
                else
                    # False positive with probability false_positive_rate
                    rand() < params.false_positive_rate
                end

                # If detected in any individual, mark the strain as detected
                if detected
                    strain_detected[strain] = true
                end
            end
        end

        # Record which strains were detected at this timestep
        detect_matrix[t, :] = strain_detected
    end

    return detect_matrix
end

"""
    virtual_ecologist_sample(;
        virtual_population,
        proportion_sampled::Float64,
        false_positive_rate::Float64,
        false_negative_rate::Float64
    ) -> Matrix{Bool}

Legacy interface for virtual ecologist sampling. Converts parameters to the new type
system and calls the `sample` function.

# Arguments
- `virtual_population`: Population data across time steps (Vector{Vector{Matrix{Bool}}} or similar).
- `proportion_sampled::Float64`: Proportion of the population sampled at each time step (0-1).
- `false_positive_rate::Float64`: Probability of detecting a strain when it's not present (0-1).
- `false_negative_rate::Float64`: Probability of not detecting a strain when it is present (0-1).

# Returns
- `Matrix{Bool}`: Detection matrix with dimensions (n_timesteps, n_strains).
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
    @assert 0 ≤ proportion_sampled ≤ 1 "Proportion sampled must be between 0 and 1."
    @assert 0 ≤ false_positive_rate ≤ 1 "False positive rate must be between 0 and 1."
    @assert 0 ≤ false_negative_rate ≤ 1 "False negative rate must be between 0 and 1."

    # Create sampling parameters
    params = SamplingParameters(
        proportion_sampled,
        false_positive_rate,
        false_negative_rate
    )

    # Convert legacy population format to new Population type
    populations = Vector{Population}(undef, length(virtual_population))
    for t in 1:length(virtual_population)
        timestep_pop = virtual_population[t]
        individuals = Vector{Individual}(undef, length(timestep_pop))

        for i in 1:length(timestep_pop)
            # We don't have ages in the legacy format, so use 0 as a default
            individuals[i] = Individual(timestep_pop[i], 0)
        end

        populations[t] = Population(individuals)
    end

    # Call the new sample_populations function
    return sample_populations(populations, params)
end
