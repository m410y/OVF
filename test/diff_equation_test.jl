using Test

include("../src/diff_equation.jl")

@testset "euler method" begin
    @test euler_method((x, y) -> 1.0, 0:0.25:1, 0.0) == [0.0, 0.25, 0.5, 0.75, 1.0]
    @test euler_method((x, y) -> x, 0:1e-2:1, 0.0)[end] ≈ 0.5 atol=1e-2
    @test euler_method((x, y) -> y, 0:1e-2:1, 1.0)[end] ≈ exp(1) atol=3e-2
    @test euler_method((x, y) -> cos(y), 0:1e-2:1, 0.0)[end] ≈ 2atan(tanh(0.5)) atol=1e-2
    @test euler_method((x, y) -> [-y[2], y[1]], 0:1e-2:1, [1.0, 0.0])[end, :] ≈ [cos(1), sin(1)] atol=1e-2
end

@testset "rk2 method 1/2" begin
    @test rk2_method((x, y) -> 1.0, 0:0.25:1, 0.0, α=1/2) == [0.0, 0.25, 0.5, 0.75, 1.0]
    @test rk2_method((x, y) -> x, 0:1e-2:1, 0.0, α=1/2)[end] ≈ 0.5 atol=1e-4
    @test rk2_method((x, y) -> y, 0:1e-2:1, 1.0, α=1/2)[end] ≈ exp(1) atol=1e-4
    @test rk2_method((x, y) -> cos(y), 0:1e-2:1, 0.0, α=1/2)[end] ≈ 2atan(tanh(0.5)) atol=1e-4
    @test rk2_method((x, y) -> [-y[2], y[1]], 0:1e-2:1, [1.0, 0.0])[end, :] ≈ [cos(1), sin(1)] atol=1e-4
end

@testset "rk2 method 1" begin
    @test rk2_method((x, y) -> 1.0, 0:0.25:1, 0.0, α=1) == [0.0, 0.25, 0.5, 0.75, 1.0]
    @test rk2_method((x, y) -> x, 0:1e-2:1, 0.0, α=1)[end] ≈ 0.5 atol=1e-4
    @test rk2_method((x, y) -> y, 0:1e-2:1, 1.0, α=1)[end] ≈ exp(1) atol=1e-4
    @test rk2_method((x, y) -> cos(y), 0:1e-2:1, 0.0, α=1)[end] ≈ 2atan(tanh(0.5)) atol=1e-4
    @test rk2_method((x, y) -> [-y[2], y[1]], 0:1e-2:1, [1.0, 0.0], α=1)[end, :] ≈ [cos(1), sin(1)] atol=1e-4
end

@testset "rk4 method" begin
    @test rk4_method((x, y) -> 1.0, 0:0.25:1, 0.0) == [0.0, 0.25, 0.5, 0.75, 1.0]
    @test rk4_method((x, y) -> x, 0:1e-2:1, 0.0)[end] ≈ 0.5 atol=1e-6
    @test rk4_method((x, y) -> y, 0:1e-2:1, 1.0)[end] ≈ exp(1) atol=1e-6
    @test rk4_method((x, y) -> cos(y), 0:1e-2:1, 0.0)[end] ≈ 2atan(tanh(0.5)) atol=1e-6
    @test rk4_method((x, y) -> [-y[2], y[1]], 0:1e-2:1, [1.0, 0.0])[end, :] ≈ [cos(1), sin(1)] atol=1e-6
end

@testset "euler implicit method" begin
    @test euler_implicit_method((x, y) -> 1, (x, y) -> 0, 0:0.25:1, 0.0) == [0.0, 0.25, 0.5, 0.75, 1.0]
    @test euler_implicit_method((x, y) -> x, (x, y) -> 0, 0:1e-2:1, 0.0)[end] ≈ 0.5 atol=1e-2
    @test euler_implicit_method((x, y) -> y, (x, y) -> 1, 0:1e-2:1, 1.0)[end] ≈ exp(1) atol=3e-2
    @test euler_implicit_method((x, y) -> cos(y), (x, y) -> -sin(y), 0:1e-2:1, 0.0)[end] ≈ 2atan(tanh(0.5)) atol=1e-2
    @test euler_implicit_method((x, y) -> [-y[2], y[1]], (x, y) -> [0 -1; 1 0], 0:1e-2:1, [1.0, 0.0])[end, :] ≈ [cos(1), sin(1)] atol=1e-2
end

@testset "poisson 1D bounded" begin
    @test poisson_1D(zero, 0:1e-2:1, (1, 0, 0), (1, 0, 1)) ≈ 0:1e-2:1 atol=1e-14
    @test poisson_1D(sin, 0:1e-2:1, (1, 0, 0), (1, 0, -sin(1))) ≈ -sin.(0:1e-2:1) atol=1e-2
end

@testset "conductivity 1D bounded" begin
    @test conductivity_1D((x, t) -> 0.0, 0:1e-1:1, 0:1e-1:1, (1, 0, 0), (1, 0, 0), zeros(11))[end, :] ≈ zeros(11) atol=1e-14
    @test conductivity_1D((x, t) -> 0.0, 0:1e-2:1, 0:1e-2:0.1, (1, 0, 0), (1, 0, 0), sinpi.(0:1e-2:1))[end, :] ≈ exp.(-0.1*pi^2)*sinpi.(0:1e-2:1) atol=1e-2
end