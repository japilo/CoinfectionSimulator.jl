using LinearAlgebra

@testset "Data Preparation Tests" begin
    @testset "Create Interaction Matrix" begin
        # Test basic functionality
        Random.seed!(123)

        # Test with priority effects (asymmetric matrix)
        matrix1 = create_interaction_matrix(3, true, 0.3, cf_ratio=0.5)
        @test size(matrix1) == (3, 3)
        @test all(diag(matrix1) .== 1.0)  # Diagonal should be 1.0
        @test !issymmetric(matrix1)  # Should be asymmetric with priority effects

        # Test without priority effects (symmetric matrix)
        matrix2 = create_interaction_matrix(4, false, 0.2, cf_ratio=0.6)
        @test size(matrix2) == (4, 4)
        @test all(diag(matrix2) .== 1.0)
        @test issymmetric(matrix2)  # Should be symmetric without priority effects

        # Test with zero interaction strength (should be identity matrix)
        matrix3 = create_interaction_matrix(2, true, 0.0)
        @test matrix3 ≈ Matrix{Float64}(I, 2, 2)

        # Test with high facilitation ratio
        matrix4 = create_interaction_matrix(5, false, 0.5, cf_ratio=0.9)
        off_diag = matrix4[.!LinearAlgebra.I(5)]
        @test count(x -> x > 1.0, off_diag) > count(x -> x < 1.0, off_diag)

        # Test with high competition ratio
        matrix5 = create_interaction_matrix(5, false, 0.5, cf_ratio=0.1)
        off_diag = matrix5[.!LinearAlgebra.I(5)]
        @test count(x -> x < 1.0, off_diag) > count(x -> x > 1.0, off_diag)
    end

    @testset "Parameter Validation" begin
        # Test invalid parameters
        @test_throws ArgumentError create_interaction_matrix(0, true, 0.5)  # Invalid strain count
        @test_throws ArgumentError create_interaction_matrix(3, true, -0.1)  # Negative interaction strength
        @test_throws ArgumentError create_interaction_matrix(3, true, 1.1)   # Interaction strength > 1
        @test_throws ArgumentError create_interaction_matrix(3, true, 0.5, cf_ratio=-0.1)  # Negative cf_ratio
        @test_throws ArgumentError create_interaction_matrix(3, true, 0.5, cf_ratio=1.1)   # cf_ratio > 1
    end

    @testset "DataFrame Interface" begin
        # Test DataFrame interface for creating multiple matrices
        df = DataFrame(
            interaction_strength=[0.1, 0.2, 0.3],
            cf_ratio=[0.3, 0.5, 0.7],
            priority_effects=[true, false, true],
            strains=[2, 3, 4]
        )

        matrices = create_interaction_matrix(df)

        @test length(matrices) == 3
        @test size(matrices[1]) == (2, 2)
        @test size(matrices[2]) == (3, 3)
        @test size(matrices[3]) == (4, 4)

        # Test with missing columns
        df_missing = DataFrame(
            interaction_strength=[0.1, 0.2],
            cf_ratio=[0.3, 0.5],
            # Missing priority_effects
            strains=[2, 3]
        )

        @test_throws ArgumentError create_interaction_matrix(df_missing)
    end

    @testset "Backward Compatibility" begin
        # Test that the legacy interface still works

        df = DataFrame(
            interaction_strength=[0.2, 0.0],
            cf_ratio=[0.4, 0.5],
            priority_effects=[true, false],
            strains=[3, 2]
        )

        matrices_new = create_interaction_matrix(df)
        matrices_old = prep_interaction_matrix(df)

        # Test that both functions produce same number of matrices
        @test length(matrices_new) == length(matrices_old)
        @test length(matrices_new) == 2

        # Test that dimensions match
        @test size(matrices_new[1]) == size(matrices_old[1]) == (3, 3)
        @test size(matrices_new[2]) == size(matrices_old[2]) == (2, 2)

        # Test that zero interaction matrices are identical (deterministic)
        @test matrices_new[2] == matrices_old[2]
        @test matrices_new[2] ≈ Matrix{Float64}(I, 2, 2)

        # Test single matrix creation with zero interaction (deterministic)
        matrix_new = create_interaction_matrix(3, true, 0.0, cf_ratio=0.5)
        matrix_old = prep_interaction_matrix(3, true, 0.0, cf_ratio=0.5)

        @test matrix_new == matrix_old
        @test matrix_new ≈ Matrix{Float64}(I, 3, 3)
    end
end
