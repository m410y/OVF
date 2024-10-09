using Test

include("../src/diff_equation.jl")

@testset "poisson 1D bounded" begin
    @test poisson_1D(zero, 0:1e-2:1, (1, 0, 0), (1, 0, 1)) ≈ 0:1e-2:1 atol = 1e-14
    @test poisson_1D(sin, 0:1e-2:1, (1, 0, 0), (1, 0, -sin(1))) ≈ -sin.(0:1e-2:1) atol =
        1e-2
end

@testset "conductivity 1D bounded" begin
    @test conductivity_1D(
        (x, t) -> 0.0,
        0:1e-1:1,
        0:1e-1:1,
        (1, 0, 0),
        (1, 0, 0),
        zeros(11),
    )[
        end,
        :,
    ] ≈ zeros(11) atol = 1e-14
    @test conductivity_1D(
        (x, t) -> 0.0,
        0:1e-2:1,
        0:1e-2:0.1,
        (1, 0, 0),
        (1, 0, 0),
        sinpi.(0:1e-2:1),
    )[
        end,
        :,
    ] ≈ exp.(-0.1 * pi^2) * sinpi.(0:1e-2:1) atol = 1e-2
end


@testset "schrödinger stationary 1D" begin
    @test schrödinger_stationary_1D(x -> 0, -1:1e-3:1)[1] ≈ pi^2/8 atol=1e-2
    @test schrödinger_stationary_1D(x -> x^2/2, -6:1e-2:6)[1] ≈ 0.5 atol=1e-4
end