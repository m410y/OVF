using Test

include("../src/specials.jl")

@testset "erf" begin
    @test erf(0.0) == 0.0
    @test erf(99.0) == 1.0
    @test erf(1.0) ≈ 0.8427007929497149 atol=1e-14
end

@testset "j0" begin
    @test j0(0.0) == 1.0
    @test j0(0.1) ≈ 0.99750156206604 atol=1e-14
    @test j0(1.0) ≈ 0.7651976865579666 atol=1e-14
    @test j0(2.404825557695773) ≈ 0.0 atol=1e-14
end

@testset "j1" begin
    @test j1(0.0) ≈ 0.0 atol=1e-14
    @test j1(0.1) ≈ 0.049937526036242 atol=1e-14
    @test j1(1.0) ≈ 0.44005058574493355 atol=1e-14
    @test j1(3.831705970207512) ≈ 0.0 atol=1e-14
end