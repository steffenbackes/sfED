@testset "ModelParameters" begin
	mp = sfED.ModelParameters(norb=2)
    @test typeof(mp.norb) === UInt64
    @test typeof(mp.Nstates) === UInt64
    @test typeof(mp.Nmax) === UInt64
    @test mp.norb == UInt64(2)
    @test mp.Nmax == UInt64(2*2)
    @test mp.Nstates == UInt64(4^2)
end
