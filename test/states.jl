using SparseArrays

@testset "basic properties" begin
    @test_throws DomainError sfED._nnmax(-1,1)
    @test_throws DomainError sfED._nnmax(1,-1)
    @test_throws DomainError sfED._nnmax(3,1)

    @test typeof(sfED.noSpinConfig(1,1)) === Int64   
    @test sfED.noSpinConfig(1,1) == 1
    @test_throws DomainError sfED.noSpinConfig(-1,1)
    @test_throws DomainError sfED.noSpinConfig(1,-1)
    @test_throws DomainError sfED.noSpinConfig(3,1)


    @test_throws DomainError sfED.spinConfig(0,1,1)
    @test_throws DomainError sfED.spinConfig(1,-1,1)
    @test_throws DomainError sfED.spinConfig(1,1,-1)
    @test_throws DomainError sfED.spinConfig(3,1,1)
    @test sfED.spinConfig(1,2,4) == -2
    @test sfED.spinConfig(2,2,4) == 0
    @test sfED.spinConfig(3,2,4) == 2

    @test sfED.indexSpinConfig(sfED.spinConfig(1,2,4),2,4) == 1
    @test sfED.indexSpinConfig(sfED.spinConfig(2,2,4),2,4) == 2
    @test sfED.indexSpinConfig(sfED.spinConfig(3,2,4),2,4) == 3

    @test sfED.getSpin(Int8[1]) == 1
    @test sfED.getSpin(Int8[0,1]) == -1
    @test sfED.getSpin(Int8[1,1,1,0,1]) == 2

    @test sfED.getCsign(1,Int8[1, 1, 1]) == 1
    @test sfED.getCsign(2,Int8[1, 1, 1]) == -1
    @test sfED.getCsign(3,Int8[1, 1, 1]) == 1
end
