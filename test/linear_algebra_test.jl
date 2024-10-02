using Test

include("../src/linear_algebra.jl")

@testset "identity matrix" begin
    @test identity_matrix(1) == 1
    @test identity_matrix(2) == [1 0; 0 1]
end

@testset "thomas non-volatile" begin
    @test thomas_algorithm([0., 0.], [1., 1., 1.], [0., 0.], [1., 1., 1.]) == [1., 1., 1.]
    @test thomas_algorithm([1., 1., 0.], [1., -2., -2., 1.], [0., 1., 1.], [1., 1., 1., 1.]) == [1., 0., 0., 1.]
end

@testset "thomas volatile" begin
    @test thomas_algorithm!([0., 0.], [1., 1., 1.], [0., 0.], [1., 1., 1.]) == [1., 1., 1.]
    @test thomas_algorithm!([1., 1., 0.], [1., -2., -2., 1.], [0., 1., 1.], [1., 1., 1., 1.]) == [1., 0., 0., 1.]
end