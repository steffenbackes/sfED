@testset "ModelParameters" begin
   mp = sfED.ModelParameters(norb=2)
   @test typeof(mp.norb) === Int64
   @test typeof(mp.Nstates) === Int64
   @test typeof(mp.Nmax) === Int64
   @test mp.norb == Int64(2)
   @test mp.Nmax == Int64(2*2)
   @test mp.Nstates == Int64(4^2)
end
