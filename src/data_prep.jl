"""
Data preparation utilities for coinfection simulation experiments.
"""

using Random
using Distributions
using DataFrames

"""
    prep_interaction_matrix(; df::DataFrame) -> Vector{Matrix{Float64}}

Generates interaction matrices for coinfection simulations based on DataFrame parameters.

# Arguments
- `df::DataFrame`: DataFrame containing the following required columns:
  - `:interaction_strength` (`Float64`): Standard deviation for interaction value generation
  - `:cf_ratio` (`Float64`): Mean ratio of competitive to facilitative interactions
  - `:priority_effects` (`Bool`): If true, creates asymmetric matrices; if false, symmetric
  - `:strains` (`Int`): Number of strains (determines matrix dimensions)

# Returns
- `Vector{Matrix{Float64}}`: Vector of interaction matrices, one per DataFrame row.
  Each matrix is `strains × strains` with diagonal elements = 1.0.

# Details
- Diagonal elements are always 1.0 (self-interaction)
- For `priority_effects = true`: Off-diagonal elements are independently sampled
- For `priority_effects = false`: Matrix is symmetric (M[i,j] = M[j,i])
- Values are sampled from Normal(cf_ratio, interaction_strength)

# Examples
```julia
using DataFrames

# Create parameter DataFrame
df = DataFrame(
    interaction_strength = [0.1, 0.2],
    cf_ratio = [0.8, 1.2],
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
        # Validate strain count
        if row.strains <= 0
            error("Number of strains must be positive, got: $(row.strains)")
        end
        
        # Initialize matrix with ones on diagonal
        int_matrix = zeros(Float64, row.strains, row.strains)
        for k in 1:row.strains
            int_matrix[k, k] = 1.0
        end
        
        # Fill off-diagonal elements
        for i in 1:row.strains
            for j in 1:row.strains
                if i != j  # Skip diagonal elements
                    if row.priority_effects
                        # Asymmetric: each element independently sampled
                        int_matrix[i, j] = rand(Normal(row.cf_ratio, row.interaction_strength))
                    else
                        # Symmetric: only fill upper triangle, then mirror
                        if i < j
                            value = rand(Normal(row.cf_ratio, row.interaction_strength))
                            int_matrix[i, j] = value
                            int_matrix[j, i] = value
                        end
                    end
                end
            end
        end

        int_matrix_list[row_idx] = int_matrix
    end

    return int_matrix_list
end

"""
    prep_interaction_matrix(strains::Int, matrix_type::String, interaction_strength::Float64; cf_ratio::Float64=1.0)

Generates a single interaction matrix with specified parameters (convenience function).

# Arguments
- `strains::Int`: Number of strains (matrix size will be strains × strains)
- `matrix_type::String`: Either "symmetric" or "asymmetric" 
- `interaction_strength::Float64`: Standard deviation for random interaction values
- `cf_ratio::Float64`: Mean for interaction values (default: 1.0)

# Returns
- `Matrix{Float64}`: An interaction matrix of size strains × strains

# Examples
```julia
# Create symmetric 3x3 matrix
matrix = prep_interaction_matrix(3, "symmetric", 0.5)

# Create asymmetric 2x2 matrix with custom mean
matrix = prep_interaction_matrix(2, "asymmetric", 0.3, cf_ratio=1.2)
```
"""
function prep_interaction_matrix(
    strains::Int, 
    matrix_type::String, 
    interaction_strength::Float64; 
    cf_ratio::Float64=1.0
)
    @assert strains > 0 "Number of strains must be positive"
    @assert matrix_type in ["symmetric", "asymmetric"] "Matrix type must be 'symmetric' or 'asymmetric'"
    
    # Create the matrix
        int_matrix = zeros(Float64, row.strains, row.strains)
        for k in 1:row.strains
            int_matrix[k, k] = 1.0
        end
    
    if matrix_type == "symmetric"
        # Fill upper triangle, then mirror to lower triangle
        for i in 1:strains
            for j in i+1:strains
                value = rand(Normal(cf_ratio, interaction_strength))
                int_matrix[i, j] = value
                int_matrix[j, i] = value  # Mirror for symmetry
            end
        end
    else  # asymmetric
        # Fill all off-diagonal elements independently
        for i in 1:strains
            for j in 1:strains
                if i != j
                    int_matrix[i, j] = rand(Normal(cf_ratio, interaction_strength))
                end
            end
        end
    end
    
    return int_matrix
end
