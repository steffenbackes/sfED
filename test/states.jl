include("../src/params.jl")
include("../src/states.jl")

@testset "basic properties" begin
    @test typeof(noSpinConfig(1,1)) === Int64   
    @test noSpinConfig(1,1) == 1   
end
