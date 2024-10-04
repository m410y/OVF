using Test

include("../src/derivation.jl")

@testset "derivative" begin
    @test derivative(ones(101), 1e-2) == zeros(101)
    @test derivative(0:1e-2:1.0, 1e-2) ≈ ones(101) atol = 1e-15
    @test derivative(sin.(0:1e-2:1.0), 1e-2) ≈ cos.(0:1e-2:1.0) atol = 1e-3
end

@testset "divided difference" begin
    @test divided_diff([1, 1], [0, 1]) == 0
    @test divided_diff([0, 1], [0, 1]) == 1
    @test divided_diff([0, 1, 4], [0, 1, 2]) == 1
end
