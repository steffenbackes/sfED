import sfED
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

@testset "IO" begin
   include("IO.jl")
end

@testset "example_calculations" begin
   include("example4orbDimer.jl")
   include("example6orbTrimer.jl")
end
