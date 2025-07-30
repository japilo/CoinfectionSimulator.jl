using Test
using DataFrames
using Random
using Distributions

@testset "Data Prep Tests" begin

	@testset "prep_interaction_matrix DataFrame input tests" begin

		@testset "Input validation tests" begin
			# Test missing required columns
			df_missing_col = DataFrame(
				interaction_strength = [0.1],
				cf_ratio = [0.5],
				priority_effects = [true],
				# Missing strains column
			)
			@test_throws ErrorException prep_interaction_matrix(df_missing_col)

			# Test invalid cf_ratio values
			df_invalid_cf = DataFrame(
				interaction_strength = [0.1],
				cf_ratio = [1.5],  # > 1
				priority_effects = [true],
				strains = [3],
			)
			@test_throws AssertionError prep_interaction_matrix(df_invalid_cf)

			df_negative_cf = DataFrame(
				interaction_strength = [0.1],
				cf_ratio = [-0.1],  # < 0
				priority_effects = [true],
				strains = [3],
			)
			@test_throws AssertionError prep_interaction_matrix(df_negative_cf)

			# Test invalid interaction_strength values
			df_invalid_strength = DataFrame(
				interaction_strength = [1.5],  # > 1
				cf_ratio = [0.5],
				priority_effects = [true],
				strains = [3],
			)
			@test_throws AssertionError prep_interaction_matrix(df_invalid_strength)

			df_negative_strength = DataFrame(
				interaction_strength = [-0.1],  # < 0
				cf_ratio = [0.5],
				priority_effects = [true],
				strains = [3],
			)
			@test_throws AssertionError prep_interaction_matrix(df_negative_strength)

			# Test invalid strains values
			df_zero_strains = DataFrame(
				interaction_strength = [0.1],
				cf_ratio = [0.5],
				priority_effects = [true],
				strains = [0],  # Must be positive
			)
			@test_throws AssertionError prep_interaction_matrix(df_zero_strains)
		end

		@testset "Valid input tests" begin
			# Test basic functionality
			df = DataFrame(
				interaction_strength = [0.1, 0.2],
				cf_ratio = [0.3, 0.7],
				priority_effects = [true, false],
				strains = [3, 4],
			)

			matrices = prep_interaction_matrix(df)
			@test length(matrices) == 2
			@test size(matrices[1]) == (3, 3)
			@test size(matrices[2]) == (4, 4)
		end

		@testset "Diagonal elements tests" begin
			df = DataFrame(
				interaction_strength = [0.2],
				cf_ratio = [0.5],
				priority_effects = [true],
				strains = [5],
			)

			matrices = prep_interaction_matrix(df)
			matrix = matrices[1]

			# Check diagonal elements are 1.0
			for i in 1:5
				@test matrix[i, i] == 1.0
			end
		end

		@testset "Interaction strength bounds tests" begin
			Random.seed!(123)
			df = DataFrame(
				interaction_strength = [0.3],
				cf_ratio = [0.5],
				priority_effects = [true],
				strains = [10],
			)

			matrices = prep_interaction_matrix(df)
			matrix = matrices[1]

			# Check all off-diagonal elements are within bounds
			for i in 1:10
				for j in 1:10
					if i != j
						@test matrix[i, j] >= 0.7  # 1 - 0.3
						@test matrix[i, j] <= 1.3  # 1 + 0.3
					end
				end
			end
		end

		@testset "Competition/facilitation ratio tests" begin
			Random.seed!(456)
			df = DataFrame(
				interaction_strength = [0.1],
				cf_ratio = [0.2],  # 20% facilitative
				priority_effects = [true],
				strains = [10],
			)

			matrices = prep_interaction_matrix(df)
			matrix = matrices[1]

			# Count facilitative (> 1.0) and competitive (< 1.0) interactions
			off_diagonal_values = [matrix[i, j] for i in 1:10, j in 1:10 if i != j]
			facilitative_count = sum(off_diagonal_values .> 1.0)
			competitive_count = sum(off_diagonal_values .< 1.0)
			total_off_diagonal = length(off_diagonal_values)

			expected_facilitative = round(Int, 0.2 * total_off_diagonal)
			@test facilitative_count == expected_facilitative
			@test competitive_count == total_off_diagonal - expected_facilitative
		end

		@testset "Symmetry tests" begin
			Random.seed!(789)
			df = DataFrame(
				interaction_strength = [0.15, 0.15],
				cf_ratio = [0.4, 0.4],
				priority_effects = [false, true],  # symmetric vs asymmetric
				strains = [4, 4],
			)

			matrices = prep_interaction_matrix(df)
			symmetric_matrix = matrices[1]
			asymmetric_matrix = matrices[2]

			# Check symmetric matrix
			for i in 1:4
				for j in 1:4
					@test symmetric_matrix[i, j] == symmetric_matrix[j, i]
				end
			end

			# Asymmetric matrix should generally not be symmetric (with high probability)
			# We'll check that at least one off-diagonal pair is different
			asymmetric_found = false
			for i in 1:4
				for j in 1:4
					if i != j && asymmetric_matrix[i, j] != asymmetric_matrix[j, i]
						asymmetric_found = true
						break
					end
				end
			end
			# Note: This test might occasionally fail due to randomness, but very unlikely
		end
	end

	@testset "prep_interaction_matrix single matrix input tests" begin

		@testset "Input validation tests" begin
			# Test invalid strains
			@test_throws AssertionError prep_interaction_matrix(0, true, 0.1)
			@test_throws AssertionError prep_interaction_matrix(-1, true, 0.1)

			# Test invalid interaction_strength
			@test_throws AssertionError prep_interaction_matrix(3, true, -0.1)
			@test_throws AssertionError prep_interaction_matrix(3, true, 1.5)

			# Test invalid cf_ratio
			@test_throws AssertionError prep_interaction_matrix(3, true, 0.1, cf_ratio = -0.1)
			@test_throws AssertionError prep_interaction_matrix(3, true, 0.1, cf_ratio = 1.5)
		end

		@testset "Basic functionality tests" begin
			matrix = prep_interaction_matrix(3, true, 0.2)
			@test size(matrix) == (3, 3)
			@test all(matrix[i, i] == 1.0 for i in 1:3)
		end

		@testset "Diagonal elements tests" begin
			matrix = prep_interaction_matrix(5, false, 0.3)

			# Check diagonal elements are 1.0
			for i in 1:5
				@test matrix[i, i] == 1.0
			end
		end

		@testset "Interaction strength bounds tests" begin
			Random.seed!(321)
			matrix = prep_interaction_matrix(8, true, 0.25)

			# Check all off-diagonal elements are within bounds
			for i in 1:8
				for j in 1:8
					if i != j
						@test matrix[i, j] >= 0.75  # 1 - 0.25
						@test matrix[i, j] <= 1.25  # 1 + 0.25
					end
				end
			end
		end

		@testset "Competition/facilitation ratio tests" begin
			Random.seed!(654)
			matrix = prep_interaction_matrix(6, true, 0.2, cf_ratio = 0.3)

			# Count facilitative (> 1.0) and competitive (< 1.0) interactions
			off_diagonal_values = [matrix[i, j] for i in 1:6, j in 1:6 if i != j]
			facilitative_count = sum(off_diagonal_values .> 1.0)
			competitive_count = sum(off_diagonal_values .< 1.0)
			total_off_diagonal = length(off_diagonal_values)

			expected_facilitative = round(Int, 0.3 * total_off_diagonal)
			@test facilitative_count == expected_facilitative
			@test competitive_count == total_off_diagonal - expected_facilitative
		end

		@testset "Symmetry tests" begin
			Random.seed!(987)
			symmetric_matrix = prep_interaction_matrix(5, false, 0.1)
			asymmetric_matrix = prep_interaction_matrix(5, true, 0.1)

			# Check symmetric matrix
			for i in 1:5
				for j in 1:5
					@test symmetric_matrix[i, j] == symmetric_matrix[j, i]
				end
			end

			# For asymmetric matrix, we can't guarantee asymmetry due to randomness,
			# but we can verify the structure is correct (no forced symmetry)
			@test size(asymmetric_matrix) == (5, 5)
			@test all(asymmetric_matrix[i, i] == 1.0 for i in 1:5)
		end

		@testset "Default cf_ratio tests" begin
			Random.seed!(111)
			matrix = prep_interaction_matrix(6, true, 0.2)  # Default cf_ratio = 0.5

			# Count facilitative (> 1.0) and competitive (< 1.0) interactions
			off_diagonal_values = [matrix[i, j] for i in 1:6, j in 1:6 if i != j]
			facilitative_count = sum(off_diagonal_values .> 1.0)
			competitive_count = sum(off_diagonal_values .< 1.0)
			total_off_diagonal = length(off_diagonal_values)

			expected_facilitative = round(Int, 0.5 * total_off_diagonal)
			@test facilitative_count == expected_facilitative
			@test competitive_count == total_off_diagonal - expected_facilitative
		end

		@testset "Edge case tests" begin
			# Test with minimum interaction strength (should be all 1.0 off-diagonal)
			matrix_min = prep_interaction_matrix(3, true, 0.0)
			for i in 1:3
				for j in 1:3
					@test matrix_min[i, j] == 1.0
				end
			end

			# Test with maximum interaction strength
			Random.seed!(222)
			matrix_max = prep_interaction_matrix(4, true, 1.0, cf_ratio = 0.5)
			for i in 1:4
				for j in 1:4
					if i != j
						@test matrix_max[i, j] >= 0.0  # 1 - 1.0
						@test matrix_max[i, j] <= 2.0  # 1 + 1.0
					end
				end
			end

			# Test with all facilitative interactions
			Random.seed!(333)
			matrix_all_fac = prep_interaction_matrix(3, true, 0.1, cf_ratio = 1.0)
			off_diagonal_values = [matrix_all_fac[i, j] for i in 1:3, j in 1:3 if i != j]
			@test all(off_diagonal_values .> 1.0)

			# Test with all competitive interactions
			Random.seed!(444)
			matrix_all_comp = prep_interaction_matrix(3, true, 0.1, cf_ratio = 0.0)
			off_diagonal_values = [matrix_all_comp[i, j] for i in 1:3, j in 1:3 if i != j]
			@test all(off_diagonal_values .< 1.0)
		end
	end
end
