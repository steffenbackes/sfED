module sfED
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
	U = 3.0
	J = 0.3
	Up = U-2*J
	t = 1.0
	mu = (U+Up+Up-J)/2      # half filling
	beta = 40.0

	pModel = ModelParameters(norb=4)
	pSimulation = SimulationParameters(U=U,J=J,t=t,mu=mu, beta=beta)
	pFreq = FrequencyMeshes(nw=501,
	                        wmin=-8.0, wmax=8.0,
					 			   iwmax=80.0,
								   beta=pSimulation.beta)

	pNumerics = NumericalParameters(delta=0.03, cutoff=1e-6, nevalsPerSubspace=100, nevalsTotalMax=400)

	println( "We have $(pModel.norb) Orbitals, #$(pModel.Nstates) states and $(pModel.Nmax) max. number of electrons" )
	
	#######################################################################
	# Main part of the Program ############################################
	#######################################################################
	
	eps             = getEps(pNumerics,pModel)  # getEps takes a small number as argument and adds random noise to the local levels to lift degeneracies and improve numerical stability
	tmatrix         = getTmatrix(pModel,pSimulation)
	Umatrix,Jmatrix = getUJmatrix(pModel,pSimulation)
	
	println("Create noninteracting single-particle Green's function...")
	gf0_w  = getG0(eps,tmatrix,pSimulation,pFreq.wf .+ im*pNumerics.delta )      # real frequencies
	gf0_iw = getG0(eps,tmatrix,pSimulation,im*pFreq.iwf )              # Matsubara frequencies
	writeGF("gf0_w.dat",gf0_w,pFreq.wf)
	writeGF("gf0_iw.dat",gf0_iw, pFreq.iwf )
	
	allstates = generateStates(pModel)                                          # generate all Fock states as arrays
	#printStateInfo(allstates)
	
	evallist,eveclist = getEvalveclist(eps,tmatrix,Umatrix,Jmatrix,pSimulation.mu,allstates,pNumerics)   # Setup Hamiltonian and solve it
	
	println("Groundstate energy E0=", minimum( first.(evallist) )  )
	println("Partition function Z=",getZ(evallist,pSimulation.beta) )
	printEvalInfo(evallist,eveclist,allstates)
	
	println("Create interacting single-particle Green's function...")
	gf_w, gf_iw, evalContributions = getGF(evallist,eveclist,allstates,pModel,pSimulation,pFreq,pNumerics)
	
	sigma_w    = getSigma(gf0_w,gf_w)                                     # get Selfenergy
	sigma_iw   = getSigma(gf0_iw,gf_iw)
	
	writeGF("gf_w.dat",    gf_w,    pFreq.wf)                        # write output
	writeGF("gf_iw.dat",   gf_iw,   pFreq.iwf)
	writeGF("sigma_w.dat", sigma_w, pFreq.wf)
	writeGF("sigma_iw.dat",sigma_iw,pFreq.iwf)
	writeEvalContributionsSectors("evalContributionsSectors.dat", evalContributions)
	writeEvalContributions("evalContributions.dat", evalContributions)
	end # end example function
end
