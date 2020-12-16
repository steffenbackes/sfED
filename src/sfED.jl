module sfED
__precompile__(false)

export example_run, noSpinConfig

using Arpack
using LinearAlgebra
using SparseArrays
using Printf
using Random

include("params.jl")
include("states.jl")
include("hamiltonian.jl")
include("IO.jl")
include("greensfunction.jl")

function example_run()
   norb = 4
   U = 3.0
   J = 0.3
   Up = U-2*J
   t = 1.0
   mu = (U+Up+Up-J)/2      # half filling
   beta = 40.0
   gf_flav = [1,5]

   pSimulation = SimulationParameters(U=U,Up=Up,J=J,t=t,mu=mu, beta=beta, gf_flav=gf_flav)
   pFreq = FrequencyMeshes(nw=501,
                           wmin=-8.0, wmax=8.0,
                           iwmax=80.0,
                           beta=pSimulation.beta)

   pNumerics = NumericalParameters(delta=0.03, cutoff=1e-6)

   # generate all the states
   fockstates = Fockstates(norb=norb)

   println( "We have $(fockstates.norb) Orbitals, #$(fockstates.Nstates) states and $(fockstates.Nmax) max. number of electrons" )
   
   #######################################################################
   # Main part of the Program ############################################
   #######################################################################
   
   eps             = getEps(fockstates,pNumerics,)  # getEps takes a small number as argument and adds random noise to the local levels to lift degeneracies
   tmatrix         = getTmatrix(fockstates,pSimulation)
   Umatrix,Jmatrix = getUJmatrix(fockstates,pSimulation)

   println("Create noninteracting single-particle Green's function...")
   gf0_w  = getG0(eps,tmatrix,pSimulation,FrequencyMeshCplx(pFreq.wf .+ im*pNumerics.delta) )    # real frequencies
   gf0_iw = getG0(eps,tmatrix,pSimulation,FrequencyMeshCplx(im*pFreq.iwf) )                      # Matsubara frequencies
   writeGF("gf0_w.dat",gf0_w,pFreq.wf)
   writeGF("gf0_iw.dat",gf0_iw, pFreq.iwf )
   
   #writeStateInfo(allstates)
   
   eigenspace = Eigenspace(eps,tmatrix,Umatrix,Jmatrix,pSimulation.mu,fockstates,pNumerics)   # Setup Hamiltonian and solve it, result is ordered by N,S
   
   println("Groundstate energy E0=", minimum( eigenspace.evals  ) )
   println("Partition function Z=",getZ(eigenspace.evals,pSimulation.beta) )
   writeEvalInfo(eigenspace,fockstates)

   println("Determining overlaps between eigenvectors...")
   transitions1pGF = getPossibleTransitions(eigenspace,fockstates,pSimulation.gf_flav,pSimulation.beta,pNumerics,1)   # contains list of possible transitions, E1,E2 and the overlap elements
#   writeTransitionsOverlaps("transitionOverlaps.dat",overlaps1pGF) # This file gets HUUGE!!
#
   println("Create interacting single-particle Green's function...")
   gf_w, gf_iw = getGF(transitions1pGF,getZ(eigenspace.evals,pSimulation.beta),pSimulation,pFreq,pNumerics)

   sigma_w    = getSigma(gf0_w,gf_w)                                     # get Selfenergy
   sigma_iw   = getSigma(gf0_iw,gf_iw)
   
   writeGF("gf_w.dat",    gf_w,    pFreq.wf)                        # write output
   writeGF("gf_iw.dat",   gf_iw,   pFreq.iwf)
   writeGF("sigma_w.dat", sigma_w, pFreq.wf)
   writeGF("sigma_iw.dat",sigma_iw,pFreq.iwf)
#   writeEvalContributionsSectors("evalContributionsSectors.dat", evalContributions)
#   writeEvalContributions("evalContributions.dat", evalContributions)

#  println("Create interacting two-particle Green's function...")
#  gf2part,evalContributions = getGF2part(evallist,eveclist,allstates,NSperm,pModel,pSimulation,pFreq,pNumerics)
#  writeGF2part("gf2part_w1w2.dat",   gf2part,   pFreq.iwf)
#  writeEvalContributionsSectors("eval2partContributionsSectors.dat", evalContributions)
#  writeEvalContributions("eval2partContributions.dat", evalContributions)

   end # end example function
end
#example_run()
