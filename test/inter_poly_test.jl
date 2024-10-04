using Test

include("../src/inter_poly.jl")

@testset "eval polynom" begin
    @test eval_polynom([0], 0) == 0
    @test eval_polynom([1, 2, 1], 3) == 16
    @test eval_polynom([0, 0, 1, 0, 0], -1) == 1
end

@testset "root polynom" begin
    @test root_polynom([0, 0, 0]) == [0, 0, 0, 1]
    @test root_polynom([0, 1, 2]) == [0, 2, -3, 1]
end

@testset "lagrange polynom" begin
    @test lagrange_polynom([1], [0]) == [1]
    @test lagrange_polynom([1, 1], [0, 1]) == [1, 0]
    @test lagrange_polynom([0, 1, 4], [0, 1, 2]) == [0, 0, 1]
end

@testset "newton polynom" begin
    @test newton_polynom([1], [0]) == [1]
    @test newton_polynom([1, 1], [0, 1]) == [1, 0]
    @test newton_polynom([0, 1, 4], [0, 1, 2]) == [0, 0, 1]
end
