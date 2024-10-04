using Test

include("../src/integration.jl")

@testset "left integral" begin
    @test left_integral(x -> x^2, 0.0:1e-2:1.0) ≈ 1 / 3 atol = 1e-2
    @test left_integral(x -> sin(x), 0.0:1e-2:pi/2) ≈ 1.0 atol = 1e-2
end

@testset "right integral" begin
    @test right_integral(x -> x^2, 0.0:1e-2:1.0) ≈ 1 / 3 atol = 1e-2
    @test right_integral(x -> sin(x), 0.0:1e-2:pi/2) ≈ 1.0 atol = 1e-2
end

@testset "trapeziod integral" begin
    @test trapezoid_integral(x -> x^2, 0.0:1e-2:1.0) ≈ 1 / 3 atol = 1e-4
    @test trapezoid_integral(x -> sin(x), 0.0:1e-2:pi/2) ≈ 1.0 atol = 1e-3
end

@testset "middle integral" begin
    @test middle_integral(x -> x^2, 0.0:1e-2:1.0) ≈ 1 / 3 atol = 1e-4
    @test middle_integral(x -> sin(x), 0.0:1e-2:pi/2) ≈ 1.0 atol = 1e-3
end

@testset "simpson integral" begin
    @test simpson_integral(x -> x^2, 0.0:1e-2:1.0) ≈ 1 / 3 atol = 1e-4
    @test simpson_integral(x -> sin(x), 0.0:1e-2:pi/2) ≈ 1.0 atol = 1e-3
end
