using Test
using SLATEC

@testset "SLATEC.jl" begin

    @testset "BLAS" begin
        @testset "daxpy!" begin
            y = [1.0, 2.0, 3.0]
            x = [1.0, 1.0, 1.0]
            daxpy!(2.0, x, y)
            @test y == [3.0, 4.0, 5.0]

            # Zero alpha - no change
            y = [1.0, 2.0, 3.0]
            daxpy!(0.0, x, y)
            @test y == [1.0, 2.0, 3.0]
        end

        @testset "dscal!" begin
            x = [1.0, 2.0, 3.0]
            dscal!(2.0, x)
            @test x == [2.0, 4.0, 6.0]
        end

        @testset "dcopy!" begin
            x = [1.0, 2.0, 3.0]
            y = [0.0, 0.0, 0.0]
            dcopy!(x, y)
            @test y == x
        end

        @testset "dswap!" begin
            x = [1.0, 2.0, 3.0]
            y = [4.0, 5.0, 6.0]
            dswap!(x, y)
            @test x == [4.0, 5.0, 6.0]
            @test y == [1.0, 2.0, 3.0]
        end

        @testset "drotg" begin
            # 3-4-5 triangle
            r, z, c, s = drotg(3.0, 4.0)
            @test r ≈ 5.0
            @test c ≈ 0.6
            @test s ≈ 0.8

            # Zero b
            r, z, c, s = drotg(5.0, 0.0)
            @test r == 5.0
            @test c == 1.0
            @test s == 0.0
        end
    end

    @testset "MINPACK" begin
        @testset "denorm" begin
            # Basic cases
            @test denorm([3.0, 4.0]) ≈ 5.0
            @test denorm([1.0, 0.0, 0.0]) ≈ 1.0
            @test denorm([0.0, 0.0, 0.0]) == 0.0

            # Overflow protection - large values
            large = [1e300, 1e300, 1e300]
            @test isfinite(denorm(large))
            @test denorm(large) ≈ sqrt(3) * 1e300

            # Underflow protection - small values
            small = [1e-300, 1e-300, 1e-300]
            @test denorm(small) > 0.0
            @test denorm(small) ≈ sqrt(3) * 1e-300
        end
    end

end
