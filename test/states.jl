using SparseArrays

@testset "basic properties" begin
    @test typeof(sfED.noSpinConfig(1,UInt64(1))) === UInt64   
    @test sfED.noSpinConfig(1,UInt64(1)) == UInt64(1)
end
