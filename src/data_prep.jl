"""
Data preparation utilities for coinfection simulation experiments.
"""

using Random
using Distributions
using DataFrames

"""
	prep_interaction_matrix(; df::DataFrame) -> Vector{Matrix{Float64}}

Generates interaction matrices for coinfection simulations based on DataFrame
parameters. Interactions are sampled from a uniform distribution, whose bounds are
defined by the ratio of competition to facilitation in the pathogen community, and the
interaction strength.

# Arguments
- `df::DataFrame`: DataFrame containing the following required columns:
- `:interaction_strength` (`Float64`): Defines the outer bounds of the interaction
  multipliers. For example, an interaction strength of 0.1 means that the interaction
  multipliers will be sampled from a uniform distribution between 0.9 and 1.1. Interaction
  strength must be between 0 and 1. If interaction strength is 0, all interactions will be
  set to neutral (1.0) regardless of the `cf_ratio`.
- `:cf_ratio` (`Float64`): Defines the ratio of facilitation to competition in the
  interaction matrix. A ratio less than 0.5 means that the matrix will have more
  competitive interactions, while a ratio greater than 0.5 means that the matrix will have
  more facilitative interactions. The ratio must be between 0 and 1.
- `:priority_effects` (`Bool`): When `true`, the interaction matrix is asymmetric
  (priority effects). When `false`, the matrix is symmetric (no priority effects).
- `:strains` (`Int`): Number of strains in the pathogen community (determines interaction
  matrix dimensions).

# Returns
- `Vector{Matrix{Float64}}`: Vector of interaction matrices, one per DataFrame row.
  Each matrix is `strains × strains` with diagonal elements = 1.0.

# Details
- Diagonal elements are always 1.0 (self-interaction)
- For `priority_effects = true`: Off-diagonal elements are independently sampled
- For `priority_effects = false`: Matrix is symmetric (M[i,j] = M[j,i])

# Examples
```julia
using DataFrames

# Create parameter DataFrame
df = DataFrame(
	interaction_strength = [0.1, 0.2],
	cf_ratio = [0.3, 0.7],
	priority_effects = [true, false],
	strains = [3, 4]
)

# Generate interaction matrices
matrices = prep_interaction_matrix(df)
```
"""
function prep_interaction_matrix(df::DataFrame)
	# Validate required columns
	required_columns = [:interaction_strength, :cf_ratio, :priority_effects, :strains]
	for col in required_columns
		if !hasproperty(df, col)
			error("DataFrame must contain column: $col")
		end
	end

	@assert all(cf_ratio -> cf_ratio <= 1, df.cf_ratio) "cf_ratio must be less than or equal to 1"
	@assert all(cf_ratio -> cf_ratio >= 0, df.cf_ratio) "cf_ratio must be non-negative"
	@assert all(strains -> strains > 0, df.strains) "strains must be positive integers"
	@assert all(interaction_strength -> interaction_strength >= 0, df.interaction_strength) "interaction_strength must be non-negative"
	@assert all(interaction_strength -> interaction_strength <= 1, df.interaction_strength) "interaction_strength must be less than or equal to 1"
	@assert all(priority_effects -> isa(priority_effects, Bool), df.priority_effects) "priority_effects must be Boolean"

	# Ensure correct column types
	df = DataFrame(df)  # Create a copy to avoid modifying the original
	df.interaction_strength = convert(Vector{Float64}, df.interaction_strength)
	df.cf_ratio = convert(Vector{Float64}, df.cf_ratio)
	df.priority_effects = convert(Vector{Bool}, df.priority_effects)
	df.strains = convert(Vector{Int}, df.strains)

	# Initialize result vector
	n_rows = nrow(df)
	int_matrix_list = Vector{Matrix{Float64}}(undef, n_rows)

	# Generate matrices for each row
	for (row_idx, row) in enumerate(eachrow(df))

		int_matrix = ones(Float64, row.strains, row.strains)

		# Get all off-diagonal linear indices
		all_indices = 1:(row.strains^2)
		diagonal_indices = [i + (i-1)*row.strains for i in 1:row.strains]
		off_diagonal_indices = setdiff(all_indices, diagonal_indices)

		# Sample facilitative indices
		n_facilitative = round(Int, row.cf_ratio * length(off_diagonal_indices))
		fac_indices = sample(off_diagonal_indices, n_facilitative, replace = false)
		comp_indices = setdiff(off_diagonal_indices, fac_indices)

		# Fill matrix
		if row.interaction_strength != 0
			int_matrix[fac_indices] = rand(Uniform(1.0, 1 + row.interaction_strength), length(fac_indices))
			int_matrix[comp_indices] = rand(Uniform(1 - row.interaction_strength, 1.0), length(comp_indices))
		else
			int_matrix[fac_indices] .= 1.0
			int_matrix[comp_indices] .= 1.0
		end

		# Handle symmetry if needed
		if !row.priority_effects
			# Make symmetric by copying upper triangle to lower triangle
			for i in 1:row.strains
				for j in (i+1):row.strains
					int_matrix[j, i] = int_matrix[i, j]
				end
			end
		end

		int_matrix_list[row_idx] = int_matrix
	end

	return int_matrix_list
end

"""
	prep_interaction_matrix(strains::Int, priority_effects::Bool, interaction_strength::Float64; cf_ratio::Float64=1.0)

Generates a single interaction matrix with specified parameters. Interactions are sampled from a 
  uniform distribution, whose bounds are defined by the ratio of competition to facilitation in the pathogen community,
  and the interaction strength.

# Arguments
- `strains::Int`: Number of strains (matrix size will be strains × strains)
- `priority_effects::Bool`: When `true`, the interaction matrix is asymmetric (priority effects). When `false`, the matrix is symmetric (no priority effects).
- `interaction_strength::Float64`: Defines the outer bounds of the interaction
  multipliers. For example, an interaction strength of 0.1 means that the interaction
  multipliers will be sampled from a uniform distribution between 0.9 and 1.1. Interaction
  strength must be between 0 and 1. If interaction strength is 0, all interactions will be
  set to neutral (1.0) regardless of the `cf_ratio`.
- `cf_ratio::Float64`: Defines the ratio of facilitation to competition in the
  interaction matrix. A ratio less than 0.5 means that the matrix will have more
  competitive interactions, while a ratio greater than 0.5 means that the matrix will have
  more facilitative interactions. The ratio must be between 0 and 1. The default is 0.5.

# Returns
- `Matrix{Float64}`: An interaction matrix of size strains × strains

# Examples
```julia
# Create symmetric 3x3 matrix
matrix = prep_interaction_matrix(3, "symmetric", 0.5)

# Create asymmetric 2x2 matrix with custom competition/facilitation ratio
matrix = prep_interaction_matrix(2, "asymmetric", 0.3, cf_ratio=0.8)
```
"""
function prep_interaction_matrix(
	strains::Int,
	priority_effects::Bool,
	interaction_strength::Float64;
	cf_ratio::Float64 = 0.5,
)
	@assert strains > 0 "Number of strains must be positive"
	@assert priority_effects isa Bool "priority_effects must be Boolean"
	@assert interaction_strength >= 0 "interaction_strength must be non-negative"
	@assert interaction_strength <= 1 "interaction_strength must be less than or equal to 1"
	@assert cf_ratio >= 0 "cf_ratio must be non-negative"
	@assert cf_ratio <= 1 "cf_ratio must be less than or equal to 1"

	# Create the matrix
	int_matrix = ones(Float64, strains, strains)

	# Get all off-diagonal linear indices
	all_indices = 1:(strains^2)
	diagonal_indices = [i + (i-1)*strains for i in 1:strains]
	off_diagonal_indices = setdiff(all_indices, diagonal_indices)

	# Sample facilitative indices
	n_facilitative = round(Int, cf_ratio * length(off_diagonal_indices))
	fac_indices = sample(off_diagonal_indices, n_facilitative, replace = false)
	comp_indices = setdiff(off_diagonal_indices, fac_indices)

	# Fill matrix
	if interaction_strength == 0
		# If interaction strength is 0, all interactions are neutral
		int_matrix[fac_indices] .= 1.0
		int_matrix[comp_indices] .= 1.0
	else
		# Sample from uniform distribution based on interaction strength
		int_matrix[fac_indices] = rand(Uniform(1.0, 1 + interaction_strength), length(fac_indices))
		int_matrix[comp_indices] = rand(Uniform(1 - interaction_strength, 1.0), length(comp_indices))
	end

	# Handle symmetry if needed
	if !priority_effects
		# Make symmetric by copying upper triangle to lower triangle
		for i in 1:strains
			for j in (i+1):strains
				int_matrix[j, i] = int_matrix[i, j]
			end
		end
	end

	return int_matrix
end
