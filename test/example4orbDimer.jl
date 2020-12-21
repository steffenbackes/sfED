using Arpack
using LinearAlgebra
using SparseArrays
using Random

@testset "4-orbital dimer model" begin
   norb = 4
   U = 3.0
   J = 0.3
   Up = U-2*J
   t = 1.0
   mu = (U+Up+Up-J)/2      # half filling
   beta = 40.0
   aim = 0
   gf_flav = [1]

   pSimulation = sfED.SimulationParameters(U=U,Up=Up,J=J,t=t,mu=mu, beta=beta, aim=aim, gf_flav=gf_flav)
   pFreq = sfED.FrequencyMeshes(nw=1,
                           wmin=0.68662673, wmax=1.0,
                           iwmax=2.0,
                           beta=pSimulation.beta)

   pNumerics = sfED.NumericalParameters(delta=0.03, cutoff=1e-6)

   # generate all the states
   fockstates = sfED.Fockstates(norb=norb)
   
   #######################################################################
   # Main part of the Program ############################################
   #######################################################################
  
   eps = zeros(Float64,fockstates.norb*2)
   for i=0:fockstates.norb-1              # Just shift one orbital down, the other up by +-1
      eps[2*i+1] = 1.0*(-1)^i    # up spin
      eps[2*i+2] = eps[2*i+1]    #dn spin
      # add a very small random term to each local level to lift degeneracy and improve numerical stability
      eps[2*i+1] += rand([-1,1]) * rand(Float64) * pNumerics.cutoff
      eps[2*i+2] += rand([-1,1]) * rand(Float64) * pNumerics.cutoff
   end
   tmatrix = -[0 0 t 0;
               0 0 0 t;
               t 0 0 0;
               0 t 0 0]
   Umatrix = [U Up 0 0;   # we assume two orbitals per site
              Up U 0 0;
              0 0 U Up;
              0 0 Up U]
   Jmatrix = [0 J 0 0;
              J 0 0 0;
              0 0 0 J;
              0 0 J 0]
   gf0_w  = sfED.getG0(eps,tmatrix,pSimulation,sfED.FrequencyMeshCplx(pFreq.wf .+ im*pNumerics.delta) )    # real frequencies
   gf0_iw = sfED.getG0(eps,tmatrix,pSimulation,sfED.FrequencyMeshCplx(im*pFreq.iwf) )                      # Matsubara frequencies

   eigenspace = sfED.Eigenspace(eps,tmatrix,Umatrix,Jmatrix,pSimulation,fockstates,pNumerics)   # Setup Hamiltonian and solve it, result is ordered by N,S

   transitions1pGF = sfED.get1pGFTransitions(pSimulation.gf_flav,eigenspace,fockstates,pSimulation.beta,pNumerics)   # contains list of possible transitions, E1,E2 and the overlap elements

   gf_w, gf_iw = sfED.getGF(transitions1pGF,sfED.getZ(eigenspace,pSimulation.beta),pSimulation,pFreq,pNumerics)
   sigma_w    = sfED.getSigma(gf0_w,gf_w)                                     # get Selfenergy
   sigma_iw   = sfED.getSigma(gf0_iw,gf_iw)
  
  #################################
  @test eigenspace.E0≈-13.084236
  @test gf_w[1,1,1]≈-2.7259717-0.57595420*im
  @test gf0_w[1,1,1]≈0.29096091-im*2.5399411e-03
  @test sigma_w[1,1,1]≈3.7877913-im*4.4195615E-02
end
