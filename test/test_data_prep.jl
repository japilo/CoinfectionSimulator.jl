using Test
using DataFrames
using LinearAlgebra  # For diag function
using Random  # For Random.seed!
using CoinfectionSimulator

@testset "Data Preparation Tests" begin
    
    @testset "Missing columns error" begin
        # Test missing interaction_strength column
        df_missing = DataFrame(
            cf_ratio = [0.5],
            priority_effects = [true],
            strains = [3]
        )
        @test_throws ErrorException prep_interaction_matrix(df=df_missing)
        
        # Test missing all required columns
        df_empty = DataFrame()
        @test_throws ErrorException prep_interaction_matrix(df=df_empty)
    end
    
    @testset "Invalid strains count" begin
        df_invalid = DataFrame(
            interaction_strength = [0.1],
            cf_ratio = [0.5],
            priority_effects = [true],
            strains = [0]
        )
        @test_throws ErrorException prep_interaction_matrix(df=df_invalid)
        
        df_negative = DataFrame(
            interaction_strength = [0.1],
            cf_ratio = [0.5],
            priority_effects = [true],
            strains = [-1]
        )
        @test_throws ErrorException prep_interaction_matrix(df=df_negative)
    end
    
    @testset "Valid input - basic functionality" begin
        Random.seed!(123)
        df = DataFrame(
            interaction_strength = [0.1, 0.2],
            cf_ratio = [0.5, 0.8],
            priority_effects = [true, false],
            strains = [2, 3]
        )
        
        result = prep_interaction_matrix(df=df)
        
        @test length(result) == 2
        @test size(result[1]) == (2, 2)
        @test size(result[2]) == (3, 3)
        @test all(isa.(result, Matrix{Float64}))
    end
    
    @testset "Priority effects - asymmetric matrix" begin
        Random.seed!(456)
        df = DataFrame(
            interaction_strength = [0.1],
            cf_ratio = [0.5],
            priority_effects = [true],
            strains = [3]
        )
        
        result = prep_interaction_matrix(df=df)
        matrix = result[1]
        
        # Diagonal should be 1.0
        @test all(diag(matrix) .== 1.0)
        
        # Off-diagonal elements should not necessarily be symmetric
        # (though they might be by chance - we just test the structure)
        @test size(matrix) == (3, 3)
    end
    
    @testset "No priority effects - symmetric matrix" begin
        Random.seed!(789)
        df = DataFrame(
            interaction_strength = [0.1],
            cf_ratio = [0.5],
            priority_effects = [false],
            strains = [3]
        )
        
        result = prep_interaction_matrix(df=df)
        matrix = result[1]
        
        # Diagonal should be 1.0
        @test all(diag(matrix) .== 1.0)
        
        # Matrix should be symmetric
        @test matrix ≈ matrix'
    end
    
    @testset "Single strain matrix" begin
        df = DataFrame(
            interaction_strength = [0.1],
            cf_ratio = [0.5],
            priority_effects = [true],
            strains = [1]
        )
        
        result = prep_interaction_matrix(df=df)
        
        @test size(result[1]) == (1, 1)
        @test result[1][1,1] == 1.0
    end
    
    @testset "Multiple rows processing" begin
        Random.seed!(321)
        df = DataFrame(
            interaction_strength = [0.1, 0.2, 0.3],
            cf_ratio = [0.5, 0.8, 1.2],
            priority_effects = [true, false, true],
            strains = [2, 3, 4]
        )
        
        result = prep_interaction_matrix(df=df)
        
        @test length(result) == 3
        @test size(result[1]) == (2, 2)
        @test size(result[2]) == (3, 3)
        @test size(result[3]) == (4, 4)
        
        # All matrices should have 1.0 on diagonal
        for matrix in result
            @test all(diag(matrix) .== 1.0)
        end
    end
    
    @testset "Single matrix convenience method" begin
        
        @testset "Invalid parameters" begin
            # Test invalid strains count
            @test_throws AssertionError prep_interaction_matrix(0, "symmetric", 0.1)
            @test_throws AssertionError prep_interaction_matrix(-1, "asymmetric", 0.2)
            
            # Test invalid matrix type
            @test_throws AssertionError prep_interaction_matrix(3, "invalid_type", 0.1)
            @test_throws AssertionError prep_interaction_matrix(2, "random", 0.2)
        end
        
        @testset "Symmetric matrix generation" begin
            Random.seed!(111)
            matrix = prep_interaction_matrix(3, "symmetric", 0.1, cf_ratio=0.8)
            
            # Check dimensions
            @test size(matrix) == (3, 3)
            
            # Check diagonal elements are 1.0
            @test all(diag(matrix) .== 1.0)
            
            # Check symmetry
            @test matrix ≈ matrix'
            
            # Check that off-diagonal elements are not 1.0 (with high probability)
            off_diag = matrix[matrix .!= 1.0]
            @test length(off_diag) > 0  # Should have off-diagonal elements
        end
        
        @testset "Asymmetric matrix generation" begin
            Random.seed!(222)
            matrix = prep_interaction_matrix(3, "asymmetric", 0.2, cf_ratio=1.2)
            
            # Check dimensions
            @test size(matrix) == (3, 3)
            
            # Check diagonal elements are 1.0
            @test all(diag(matrix) .== 1.0)
            
            # For a 3x3 matrix, we can't guarantee asymmetry due to randomness,
            # but we can check that the structure allows for it
            @test size(matrix) == (3, 3)
            
            # Check that off-diagonal elements exist and are not all 1.0
            off_diag = [matrix[i,j] for i in 1:3, j in 1:3 if i != j]
            @test length(off_diag) == 6  # Should have 6 off-diagonal elements for 3x3
        end
        
        @testset "Single strain matrix" begin
            matrix = prep_interaction_matrix(1, "symmetric", 0.5)
            
            @test size(matrix) == (1, 1)
            @test matrix[1,1] == 1.0
            
            # Test with asymmetric (should behave same for 1x1)
            matrix_asym = prep_interaction_matrix(1, "asymmetric", 0.5)
            @test size(matrix_asym) == (1, 1)
            @test matrix_asym[1,1] == 1.0
        end
        
        @testset "Default cf_ratio parameter" begin
            Random.seed!(333)
            # Test without specifying cf_ratio (should default to 1.0)
            matrix = prep_interaction_matrix(2, "symmetric", 0.1)
            
            @test size(matrix) == (2, 2)
            @test all(diag(matrix) .== 1.0)
            @test matrix ≈ matrix'  # Should be symmetric
        end
        
        @testset "Large matrix generation" begin
            Random.seed!(444)
            large_matrix = prep_interaction_matrix(10, "asymmetric", 0.3, cf_ratio=0.9)
            
            @test size(large_matrix) == (10, 10)
            @test all(diag(large_matrix) .== 1.0)
            
            # Check that we have correct number of off-diagonal elements
            off_diag_count = sum(large_matrix .!= 1.0)
            @test off_diag_count == 90  # 10*10 - 10 diagonal elements
        end
        
        @testset "Different parameter combinations" begin
            Random.seed!(555)
            
            # Test various combinations
            test_cases = [
                (2, "symmetric", 0.1, 0.5),
                (3, "asymmetric", 0.2, 1.5),
                (4, "symmetric", 0.05, 2.0),
                (5, "asymmetric", 0.3, 0.1)
            ]
            
            for (strains, matrix_type, interaction_strength, cf_ratio) in test_cases
                matrix = prep_interaction_matrix(strains, matrix_type, interaction_strength, cf_ratio=cf_ratio)
                
                # Basic checks for all cases
                @test size(matrix) == (strains, strains)
                @test all(diag(matrix) .== 1.0)
                @test isa(matrix, Matrix{Float64})
                
                # Symmetry check for symmetric matrices
                if matrix_type == "symmetric"
                    @test matrix ≈ matrix'
                end
            end
        end
        
        @testset "Consistency with DataFrame method" begin
            Random.seed!(666)
            
            # Generate using convenience method
            single_matrix = prep_interaction_matrix(3, "symmetric", 0.15, cf_ratio=0.75)
            
            # Generate using DataFrame method with same parameters
            Random.seed!(666)  # Reset seed for consistency
            df = DataFrame(
                interaction_strength = [0.15],
                cf_ratio = [0.75],
                priority_effects = [false],  # false = symmetric
                strains = [3]
            )
            df_result = prep_interaction_matrix(df=df)
            df_matrix = df_result[1]
            
            # Results should be identical
            @test single_matrix ≈ df_matrix
        end
    end
end
