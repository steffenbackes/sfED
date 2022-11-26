@testset "ModelParameters" begin
   sp = ModelParameters(U=1.0,J=1.5,t=2.0,mu=2.5,beta=3.0,aim=0,gf_flav=[1,0])
   sp2 = ModelParameters(U=1.0,J=1.5,Up=0.1,t=2.0,mu=2.5,beta=3.0,aim=0,gf_flav=[1,0])
   @test sp.U ≈ 1.0
   @test sp.J ≈ 1.5 
   @test sp.Up ≈ -2.0 
   @test sp.t ≈ 2.0 
   @test sp.mu≈ 2.5
   @test sp.beta ≈ 3.0
   @test sp.gf_flav == [1,0]
   @test sp2.Up ≈ 0.1
end

@testset "NumericalParameters" begin
   np = NumericalParameters(delta=1.0,cutoff=2.0)
   @test np.delta ≈ 1.0 
   @test np.cutoff ≈ 2.0 
end

@testset "FrequencyMeshes" begin
    t = FrequencyMeshes(;nw = 3, wmin=-2, wmax=3, iwmax=4, beta=1.1)
    @test length(t.wf) == 3
    @test length(t.iwf) == 2*4+1
    @test length(t.ivf) == 2*4
end
