using Test

include("../src/root_find.jl")

@testset "dichotomy root" begin
    @test dichotomy_root(x -> x, -1.0, 1.0, error=1e-5) ≈ 0.0 atol=1e-5
    @test dichotomy_root(x -> x^2 - 1.0, 0.0, 2.0, error=1e-5) ≈ 1.0 atol=1e-5
    @test_throws ArgumentError dichotomy_root(x -> x^2 - 1.0, 2.0, 4.0)
end

@testset "iteration root" begin
    @test iteration_root(x -> x, 1.0, error=1e-5, λ=0.5) ≈ 0.0 atol=1e-5
    @test_throws OverflowError iteration_root(x -> x, 1.0, error=1e-5, λ=-0.01)
    @test iteration_root(x -> x^2 - 1.0, 2.0, error=1e-5, λ=0.1) ≈ 1.0 atol=1e-4
end

@testset "newton root" begin
    @test newton_root(x -> x, x -> 1, 1.0, error=1e-5) ≈ 0.0 atol=1e-5
    @test newton_root(x -> x^2 - 1.0, x-> 2x, 2.0, error=1e-5) ≈ 1.0 atol=1e-5
end