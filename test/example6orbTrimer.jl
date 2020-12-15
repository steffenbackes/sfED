using Arpack
using LinearAlgebra
using SparseArrays
using Random

@testset "6-orbital trimer model" begin
   norb = 6
   U = 3.0
   J = 0.0
   Up = 0
   t = 1.0
   mu = (U+Up+Up-J)/2      # half filling
   beta = 40.0
   gf_flav = [1,5]

   pModel = sfED.ModelParameters(norb=norb)
   pSimulation = sfED.SimulationParameters(U=U,Up=Up,J=J,t=t,mu=mu, beta=beta, gf_flav=gf_flav)
   pFreq = sfED.FrequencyMeshes(nw=1,
                           wmin=-1.5169661, wmax=1.0,
                           iwmax=2.0,
                           beta=pSimulation.beta)

   pNumerics = sfED.NumericalParameters(delta=0.03, cutoff=1e-6, nevalsPerSubspace=400, nevalsTotalMax=5000)
   
   #######################################################################
   # Main part of the Program ############################################
   #######################################################################
  
   eps = zeros(Float64,pModel.Nmax)
   for i=0:pModel.norb-1              # Just shift one orbital down, the other up by +-1
      eps[2*i+1] = 0.0*(-1)^i    # up spin
      eps[2*i+2] = eps[2*i+1]    #dn spin
      # add a very small random term to each local level to lift degeneracy and improve numerical stability
      eps[2*i+1] += rand([-1,1]) * rand(Float64) * pNumerics.cutoff
      eps[2*i+2] += rand([-1,1]) * rand(Float64) * pNumerics.cutoff
   end
   tmatrix = -[0 0 t 0 0 0;
               0 0 0 t 0 0;
               t 0 0 0 t 0;
               0 t 0 0 0 t;
               0 0 t 0 0 0;
               0 0 0 t 0 0]
   Umatrix = [U Up 0 0 0 0; # we assume two orbitals per site
              Up U 0 0 0 0;
              0 0 U Up 0 0;
              0 0 Up U 0 0;
              0 0 0 0 U Up;
              0 0 0 0 Up U]
   Jmatrix = [0 J 0 0 0 0;
              J 0 0 0 0 0;
              0 0 0 J 0 0;
              0 0 J 0 0 0;
              0 0 0 0 0 J;
              0 0 0 0 J 0]

   gf0_w  = sfED.getG0(eps,tmatrix,pSimulation,sfED.FrequencyMeshCplx(pFreq.wf .+ im*pNumerics.delta) )    # real frequencies
   gf0_iw = sfED.getG0(eps,tmatrix,pSimulation,sfED.FrequencyMeshCplx(im*pFreq.iwf) )                      # Matsubara frequencies
   
   allstates = sfED.generateStates(pModel)                                          # generate all Fock states as arrays
   
   evallist,eveclist = sfED.getEvalveclist(eps,tmatrix,Umatrix,Jmatrix,pSimulation.mu,allstates,pNumerics)   # Setup Hamiltonian and solve it
   
   NSperm = sfED.getNSperm(evallist)                  # get permutation which sorts the Evals for N,S,E and use this for overlap and the GF to be consistent
   overlaps1pGF,possibTransitions1pGF = sfED.getPossibleTransitions(evallist,eveclist,allstates,pSimulation.gf_flav,NSperm,pNumerics)
   
   gf_w, gf_iw, evalContributions = sfED.getGF(evallist,overlaps1pGF,possibTransitions1pGF,NSperm,pModel,pSimulation,pFreq,pNumerics)

   sigma_w    = sfED.getSigma(gf0_w,gf_w)                                     # get Selfenergy
   sigma_iw   = sfED.getSigma(gf0_iw,gf_iw)
  
  #################################
  E0 = minimum(first.(evallist))
  @test E0≈-11.970696
#  @test sfED.getZ(evallist,beta)≈3.999466 too sensitive to small numerical noise
  @test gf_w[1,1,1]≈1.4431592-im*1.0092759
  @test gf_w[1,2,1]≈3.3990114-im*1.7905190
  @test gf0_w[1,1,1]≈-14.283422-im*25.255127
  @test sigma_w[1,1,1]≈1.7646194-im*0.047359131
  @test sigma_w[1,2,1]≈-1.0361185+im*0.016052775
  @test imag(gf_iw[1,1,1])≈-0.055855915
  @test real(gf_iw[1,2,1])≈3.4700558E-01
  @test sigma_iw[1,1,1]≈1.5000614E+00-im*1.4937729E-01
  @test sigma_iw[1,2,1]≈-2.8451111E+00-im*8.4704679E-06

end
