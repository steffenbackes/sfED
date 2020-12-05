using Arpack
using LinearAlgebra
using SparseArrays
using Printf
using Random
include("states.jl")
include("hamiltonian.jl")
include("io_ops.jl")
include("greensfunction.jl")

const norb = 4
const    U = 3.0
const    J = 0.3
const   Up = (U-2*J)
const    t = 1.0
const   mu =  (U + Up + Up-J)/2.0
const beta = 40.0
const delta = 0.03
const cutoff = 1e-6
const nevalsPerSubspace = 30   # how much eigenvalues to obtain for each subspace (affects ARPACK directly)
const nevalsTotalMax = 200     # how much eigenvalues to keep in total from all subspaces (affects mostly memory and Green's function)

const nw = 502         # for writing into file and plotting
const wmin = -8.0
const wmax = 8.0

const Nmax = 2*norb
const Nstates = 4^norb
println( "We have ",norb," Orbitals, #",Nstates," states and ",Nmax," max. number of electrons" )


#######################################################################
# Main part of the Program ############################################
#######################################################################

eps             = getEps(cutoff)  # getEps takes a small number as argument and adds random noise to the local levels to lift degeneracies and improve numerical stability
tmatrix         = getTmatrix()
Umatrix,Jmatrix = getUJmatrix()

println("Create noninteracting single-particle Green's function...")
gf0_w  = getG0(eps,tmatrix, LinRange(wmin,wmax,nw) .+ im*delta )      # real frequencies
gf0_iw = getG0(eps,tmatrix, ( 2 .* collect(1:nw) .-1).*(im*pi/beta) ) # Matsubara frequencies
writeGF("gf0_w.dat",gf0_w,LinRange(wmin,wmax,nw))
writeGF("gf0_iw.dat",gf0_iw, ( 2 .* collect(1:nw) .-1).*(pi/beta) )

allstates = generateStates()                                          # generate all Fock states as arrays
#printStateInfo(allstates)

evallist,eveclist = getEvalveclist(eps,tmatrix,Umatrix,Jmatrix,allstates)   # Setup Hamiltonian and solve it

println("Groundstate energy E0=", minimum( first.(evallist) )  )
println("Partition function Z=",getZ(evallist))
printEvalInfo(evallist,eveclist,allstates)

println("Create interacting single-particle Green's function...")
gf_w, gf_iw, evalContributions = getGF(evallist,eveclist, allstates)

sigma_w    = getSigma(gf0_w,gf_w)                                     # get Selfenergy
sigma_iw   = getSigma(gf0_iw,gf_iw)

writeGF("gf_w.dat",gf_w,LinRange(wmin,wmax,nw))                        # write output
writeGF("gf_iw.dat",gf_iw, ( 2 .* collect(1:nw) .-1).*(pi/beta) )
writeGF("sigma_w.dat",sigma_w,LinRange(wmin,wmax,nw))
writeGF("sigma_iw.dat",sigma_iw, ( 2 .* collect(1:nw) .-1).*(pi/beta) )
writeEvalContributionsSectors("evalContributionsSectors.dat", evalContributions)
writeEvalContributions("evalContributions.dat", evalContributions)


