using Pkg
Pkg.activate(".")
using sfED
using TimerOutputs

to = TimerOutput()

@timeit to "example1" begin
norb = 2
U = 1.0 
J = 0.0
Up = 0
t = 1/(2*sqrt(6))
mu = (U+Up+Up-J)/2      # half filling
beta = 10.0
aim = 1
#gf_flav = [1]
gf_flav = [2*m-1 for m in 1:norb]

pSimulation = ModelParameters(U=U,Up=Up,J=J,t=t,mu=mu, beta=beta, aim=aim, gf_flav=gf_flav)
pFreq = FrequencyMeshes(nw=501,
                       wmin=-8.0, wmax=8.0,
                       iwmax=80.0,
                       beta=pSimulation.beta)

pNumerics = NumericalParameters(delta=0.03, cutoff=1e-5)

fockstates = Fockstates(norb=norb)
#writeStateInfo(fockstates)

println( "We have $(fockstates.norb) Orbitals, #$(fockstates.Nstates) states and $(fockstates.norb*2) max. number of electrons" )

#######################################################################
# Main part of the Program ############################################
#######################################################################

eps             = getEps(fockstates,pNumerics,)  # getEps takes a small number as argument and adds random noise to the local levels to lift degeneracies
tmatrix         = getTmatrix(fockstates,pSimulation)
Umatrix,Jmatrix = getUJmatrix(fockstates,pSimulation)

println("Create noninteracting single-particle Green's function...")
gf0_w  = getG0(eps,tmatrix,pSimulation,FrequencyMeshCplx(pFreq.wf .+ im*pNumerics.delta) )    # real frequencies
gf0_iw = getG0(eps,tmatrix,pSimulation,FrequencyMeshCplx(im*pFreq.iwf) )                      # Matsubara frequencies
# writeGF("gf0_w.dat",gf0_w,pFreq.wf)
# writeGF("gf0_iw.dat",gf0_iw, pFreq.iwf )
#
@timeit to "Eigenspace" eigenspace = Eigenspace(eps,tmatrix,Umatrix,Jmatrix,pSimulation,fockstates,pNumerics);   # Setup Hamiltonian and solve it, result is ordered by N,S
# println("Groundstate energy E0=", eigenspace.E0 )
# println("Partition function Z=",getZ(eigenspace,pSimulation.beta) )
# println("Total Energy E=",getE(eigenspace,pSimulation.beta) )
# writeParticleNumbers( getN(eigenspace,pSimulation.beta,fockstates) )
# writeDoubleOccupations( getNN(eigenspace,pSimulation.beta, fockstates) )
# writeEvalInfo(eigenspace,fockstates)
#
# println("Determining overlaps between eigenvectors for 1partGF...")
# transitions1pGF = get1pGFTransitions(pSimulation.gf_flav,eigenspace,fockstates,pSimulation.beta,pNumerics)   # contains list of possible transitions
# #   writeTransitionsOverlaps("transitionOverlaps.dat",overlaps1pGF) # This file gets HUUGE!!
#
# println("Create interacting single-particle Green's function...")
# gf_w, gf_iw = getGF(transitions1pGF,getZ(eigenspace,pSimulation.beta),pSimulation,pFreq,pNumerics)
#
# #sigma_w    = getSigma(gf0_w,gf_w)                                     # get Selfenergy
# sigma_iw   = getSigma(gf0_iw,gf_iw)
#
# #writeGF("gf_w.dat",    gf_w,    pFreq.wf)                        # write output
# writeGF("gf_iw.dat",   gf_iw,   pFreq.iwf)
# #writeGF("sigma_w.dat", sigma_w, pFreq.wf)
# writeGF("sigma_iw.dat",sigma_iw,pFreq.iwf)
# #writeEvalContributionsSectors("evalContributionsSectors.dat", evalContributions)
# #writeEvalContributions("evalContributions.dat", evalContributions)
#
# #   println("Determining overlaps between eigenvectors for 2partGF...")
# #   transitions2pGF = get2pGFTransitions(1,eigenspace,fockstates,pSimulation.beta,pNumerics)   # contains list of possible transitions
# #   println("Create interacting two-particle Green's function...")
# #   gf2part = getGF2part(transitions2pGF,getZ(eigenspace,pSimulation.beta),pFreq,7,pNumerics)
# #   writeGF2part("gf2part_w1w2.dat",   gf2part,   pFreq.iwf)
# #   #writeEvalContributionsSectors("eval2partContributionsSectors.dat", evalContributions)
# #   #writeEvalContributions("eval2partContributions.dat", evalContributions)
# #   #Profile.print()
#
end
