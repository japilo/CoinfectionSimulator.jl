"""
Utilities for preparing simulation data and interaction matrices.
"""

using Random
using Distributions
using DataFrames
using StatsBase: sample

"""
create_interaction_matrix(strains::Int,
  priority_effects::Bool,
  interaction_strength::Float64;
  cf_ratio::Float64=0.5) -> Matrix{Float64}

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
matrix = create_interaction_matrix(2, "asymmetric", 0.3, cf_ratio=0.8)
```
"""
function create_interaction_matrix(
    strains::Int,
    priority_effects::Bool,
    interaction_strength::Float64;
    cf_ratio::Float64=0.5
)
    strains > 0 || throw(ArgumentError("Number of strains must be positive"))
    0 ≤ interaction_strength ≤ 1 || throw(ArgumentError("Interaction strength must be between 0 and 1"))
    0 ≤ cf_ratio ≤ 1 || throw(ArgumentError("Facilitation/competition ratio must be between 0 and 1"))

    # Initialize matrix with ones on diagonal
    matrix = Matrix{Float64}(I, strains, strains)

    # If interaction strength is 0, return the identity matrix
    interaction_strength ≈ 0.0 && return matrix

    # Create off-diagonal indices
    off_diag_indices = [(i, j) for i in 1:strains for j in 1:strains if i != j]

    # Sample facilitative indices
    n_facilitative = round(Int, cf_ratio * length(off_diag_indices))
    fac_indices = StatsBase.sample(off_diag_indices, n_facilitative, replace=false)
    comp_indices = setdiff(off_diag_indices, fac_indices)

    # Fill matrix with interaction values
    for (i, j) in fac_indices
        matrix[i, j] = rand(Uniform(1.0, 1 + interaction_strength))
    end

    for (i, j) in comp_indices
        matrix[i, j] = rand(Uniform(1 - interaction_strength, 1.0))
    end

    # For symmetric matrices, ensure M[i,j] = M[j,i]
    if !priority_effects
        for i in 1:strains
            for j in (i+1):strains
                matrix[j, i] = matrix[i, j]
            end
        end
    end

    return matrix
end

"""
create_interaction_matrix(df::DataFrame) -> Vector{Matrix{Float64}}

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
matrices = create_interaction_matrix(df)
```
"""
function create_interaction_matrix(df::DataFrame)
    # Validate required columns
    required_columns = [:interaction_strength, :cf_ratio, :priority_effects, :strains]
    for col in required_columns
        if !hasproperty(df, col)
            throw(ArgumentError("DataFrame must contain column: $col"))
        end
    end

    # Create matrices for each row
    matrices = map(eachrow(df)) do row
        create_interaction_matrix(
            row.strains,
            row.priority_effects,
            row.interaction_strength;
            cf_ratio=row.cf_ratio
        )
    end

    return matrices
end

# Legacy function for backward compatibility
"""
    prep_interaction_matrix(df::DataFrame) -> Vector{Matrix{Float64}}

Legacy function for backward compatibility. Calls `create_interaction_matrix`.
"""
function prep_interaction_matrix(df::DataFrame)
    create_interaction_matrix(df)
end

"""
    prep_interaction_matrix(strains::Int, priority_effects::Bool, interaction_strength::Float64;
                           cf_ratio::Float64=0.5) -> Matrix{Float64}

Legacy function for backward compatibility. Calls `create_interaction_matrix`.
"""
function prep_interaction_matrix(
    strains::Int,
    priority_effects::Bool,
    interaction_strength::Float64;
    cf_ratio::Float64=0.5
)
    create_interaction_matrix(strains, priority_effects, interaction_strength; cf_ratio=cf_ratio)
end
