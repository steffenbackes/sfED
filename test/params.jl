@testset "ModelParameters" begin

   sp = sfED.SimulationParameters(U=1.0,J=1.5,t=2.0,mu=2.5,beta=3.0,aim=0,gf_flav=[1,0])
   sp2 = sfED.SimulationParameters(U=1.0,J=1.5,Up=0.1,t=2.0,mu=2.5,beta=3.0,aim=0,gf_flav=[1,0])
   @test sp.U ≈ 1.0
   @test sp.J ≈ 1.5 
   @test sp.Up ≈ -2.0 
   @test sp.t ≈ 2.0 
   @test sp.mu≈ 2.5
   @test sp.beta ≈ 3.0
   @test sp.gf_flav == [1,0]
   @test sp2.Up ≈ 0.1

#   @test_throws DomainError sfED.NumericalParameters(delta=1.0,cutoff=2.0)
#   @test_throws DomainError sfED.NumericalParameters(delta=1.0,cutoff=2.0)
   np = sfED.NumericalParameters(delta=1.0,cutoff=2.0)
   @test np.delta ≈ 1.0 
   @test np.cutoff ≈ 2.0 
end
