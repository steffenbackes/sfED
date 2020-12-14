using SparseArrays

@testset "basic properties" begin
    @test typeof(sfED.noSpinConfig(1,Int64(1))) === Int64   
    @test sfED.noSpinConfig(1,Int64(1)) == Int64(1)
end
