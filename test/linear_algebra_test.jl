using Test

include("../src/linear_algebra.jl")

@testset "identity matrix" begin
    @test identity_matrix(1) == 1
    @test identity_matrix(2) == [1 0; 0 1]
end

@testset "thomas non-volatile" begin
    @test thomas_algorithm([0.0, 0.0], [1.0, 1.0, 1.0], [0.0, 0.0], [1.0, 1.0, 1.0]) ==
          [1.0, 1.0, 1.0]
    @test thomas_algorithm(
        [1.0, 1.0, 0.0],
        [1.0, -2.0, -2.0, 1.0],
        [0.0, 1.0, 1.0],
        [1.0, 1.0, 1.0, 1.0],
    ) == [1.0, 0.0, 0.0, 1.0]
end

@testset "thomas volatile" begin
    @test thomas_algorithm!([0.0, 0.0], [1.0, 1.0, 1.0], [0.0, 0.0], [1.0, 1.0, 1.0]) ==
          [1.0, 1.0, 1.0]
    @test thomas_algorithm!(
        [1.0, 1.0, 0.0],
        [1.0, -2.0, -2.0, 1.0],
        [0.0, 1.0, 1.0],
        [1.0, 1.0, 1.0, 1.0],
    ) == [1.0, 0.0, 0.0, 1.0]
end

@testset "power method" begin
    @test power_method([1.0 0.0; 0.0 1.0])[1] ≈ 1.0 atol = 1e-16
    @test power_method(
        [
            7.91908 8.98168 1.11756 -3.52672
            -3.01069 -3.81374 -1.36183 0.854962
            -5.40305 -6.26107 -1.07481 2.24427
            -7.74962 -3.79237 -10.6656 9.96947
        ],
    )[1] ≈ 10.000032470011144 atol = 1e-14
end

@testset "power tridiag method" begin
    @test power_tridiag_method(zeros(2), ones(3), zeros(2))[1] ≈ 1.0 atol = 1e-16
    @test power_tridiag_method(fill(-1.0, 4), fill(2.0, 5), fill(-1.0, 4))[1] ≈
          3.732050807568877 atol = 1e-14
end

@testset "inverse tridiag method" begin
    @test inverse_tridiag_method(zeros(2), ones(3), zeros(2))[1] ≈ 1.0 atol = 1e-16
    @test inverse_tridiag_method(fill(-1.0, 4), fill(2.0, 5), fill(-1.0, 4))[1] ≈
          0.26794919243112725 atol = 1e-14
end
