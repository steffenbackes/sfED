using sfED
using Test

@testset "Example Tests" begin
    @test 1 == 1
end

@testset "states" begin
    include("states.jl")
end

@testset "params" begin
    include("params.jl") 
end

@testset "sfED_productiontests" begin
   include("sfED_productiontests.jl")
end
