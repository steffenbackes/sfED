using Test
using sfED

@testset "params" begin
    include("params.jl") 
end

@testset "states" begin
    include("states.jl")
end

@testset "hamiltonian" begin
    include("hamiltonian.jl")
end

@testset "ccdagger" begin
    include("ccdagger.jl")
end

@testset "transitions" begin
    include("transitions.jl")
end

@testset "greensfunction" begin
    include("greensfunction.jl")
end

@testset "IO" begin
   include("IO.jl")
end

@testset "example_calculations" begin
   include("example4orbDimer.jl")
   #include("example6orbTrimer.jl")
end
