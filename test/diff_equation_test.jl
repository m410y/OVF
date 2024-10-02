using Test

include("../src/diff_equation.jl")

@testset "euler method" begin
    @test euler_method((x, y) -> 1.0, 0:0.25:1, 0.0) == [0.0, 0.25, 0.5, 0.75, 1.0]
    @test euler_method((x, y) -> x, 0:1e-2:1, 0.0)[end] ≈ 0.5 atol=1e-2
    @test euler_method((x, y) -> y, 0:1e-2:1, 1.0)[end] ≈ exp(1) atol=3e-2
    @test euler_method((x, y) -> cos(y), 0:1e-2:1, 0.0)[end] ≈ 2atan(tanh(0.5)) atol=1e-2
end

@testset "rk2 method 1/2" begin
    @test rk2_method((x, y) -> 1.0, 0:0.25:1, 0.0, α=1/2) == [0.0, 0.25, 0.5, 0.75, 1.0]
    @test rk2_method((x, y) -> x, 0:1e-2:1, 0.0, α=1/2)[end] ≈ 0.5 atol=1e-4
    @test rk2_method((x, y) -> y, 0:1e-2:1, 1.0, α=1/2)[end] ≈ exp(1) atol=1e-4
    @test rk2_method((x, y) -> cos(y), 0:1e-2:1, 0.0, α=1/2)[end] ≈ 2atan(tanh(0.5)) atol=1e-4
end

@testset "rk2 method 1" begin
    @test rk2_method((x, y) -> 1.0, 0:0.25:1, 0.0, α=1) == [0.0, 0.25, 0.5, 0.75, 1.0]
    @test rk2_method((x, y) -> x, 0:1e-2:1, 0.0, α=1)[end] ≈ 0.5 atol=1e-4
    @test rk2_method((x, y) -> y, 0:1e-2:1, 1.0, α=1)[end] ≈ exp(1) atol=1e-4
    @test rk2_method((x, y) -> cos(y), 0:1e-2:1, 0.0, α=1)[end] ≈ 2atan(tanh(0.5)) atol=1e-4
end

@testset "rk4 method" begin
    @test rk4_method((x, y) -> 1.0, 0:0.25:1, 0.0) == [0.0, 0.25, 0.5, 0.75, 1.0]
    @test rk4_method((x, y) -> x, 0:1e-2:1, 0.0)[end] ≈ 0.5 atol=1e-6
    @test rk4_method((x, y) -> y, 0:1e-2:1, 1.0)[end] ≈ exp(1) atol=1e-6
    @test rk4_method((x, y) -> cos(y), 0:1e-2:1, 0.0)[end] ≈ 2atan(tanh(0.5)) atol=1e-6
end