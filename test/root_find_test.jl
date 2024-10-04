using Test

include("../src/root_find.jl")

@testset "dichotomy root" begin
    @test dichotomy_root(x -> x, -1.0, 1.0) ≈ 0.0 atol = 1e-15
    @test dichotomy_root(x -> x^2 - 1.0, 0.0, 2.0) ≈ 1.0 atol = 1e-15
    @test_throws ArgumentError dichotomy_root(x -> x^2 - 1.0, 2.0, 4.0)
end

@testset "iteration root" begin
    @test iteration_root(x -> x, 1.0, λ = 0.5) ≈ 0.0 atol = 1e-15
    @test iteration_root(x -> x, 1.0, λ = -0.01, show_iter = true)[2] > 1000
    @test iteration_root(x -> x^2 - 1.0, 2.0, λ = 0.1) ≈ 1.0 atol = 1e-4
end

@testset "newton root" begin
    @test newton_root(x -> x, x -> 1, 1.0) ≈ 0.0 atol = 1e-15
    @test newton_root(x -> x^2 - 1.0, x -> 2x, 2.0) ≈ 1.0 atol = 1e-15
end

@testset "secant root" begin
    @test secant_root(x -> x, 1.0) ≈ 0.0 atol = 1e-15
    @test secant_root(x -> x^2 - 1.0, 2.0, λ = 0.1) ≈ 1.0 atol = 1e-15
end
