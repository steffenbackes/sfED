
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

    @test sfED.getSpin(sfED.FockElement[1]) == 1
    @test sfED.getSpin(sfED.FockElement[0,1]) == -1
    @test sfED.getSpin(sfED.FockElement[1,1,1,0,1]) == 2

    @test sfED.getCsign(1,sfED.FockElement[1, 1, 1]) == 1
    @test sfED.getCsign(2,sfED.FockElement[1, 1, 1]) == -1
    @test sfED.getCsign(3,sfED.FockElement[1, 1, 1]) == 1
end

@testset "generateStates" begin
    for N in 1:5
        tt = sfED.generateStates(2*N,4^N)
        # 2N+1 electrons (all sites double occupied + none occupied) for N sites
        @test length(tt) == 2N+1
        # e.g.: 1,2,3,2,1 spin configurations are possible for the 5 possible electron configurations (i.e. N = 4)
        @test all(length.(tt) .== cat(1:N,[N+1],reverse(1:N), dims=1))
    end
end
