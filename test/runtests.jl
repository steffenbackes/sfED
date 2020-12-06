using sfED
using Test

@testset "Example Tests" begin
    @test 1 == 1
end

@testset "states" begin
    include("states.jl")
end
