include("../src/params.jl")
include("../src/states.jl")
using SparseArrays

@testset "basic properties" begin
    @test typeof(noSpinConfig(1,1)) === Int64   
    @test noSpinConfig(1,1) == 1   
end
