#typedefs
const CAmatrix    = SparseMatrixCSC{Int8,Int64}   # matrix type for creation/annihilation matrix
const Hamiltonian = SparseMatrixCSC{Complex{Float32},Int64}   # matrix type for Hamiltonian matrix

const FockElement = Int8                                  # one number element of a Fock state
const Fockstate   = Array{FockElement,1}                  # The standard Fock state e.g. [00101011] 
const NSstates    = Array{Array{Array{Fockstate,1},1},1}  # An array of arrays for all Fock states for Quantum numbers N,S: typeof( allstates[N][S][i] )=Fockstate

const Eigenvalue        = Float32             
const EigenvectorElem   = Complex{Float32}
const Eigenvector       = Array{EigenvectorElem,1}
const EigenvectorMatrix = Array{Complex{Float32},2}

const FrequencyMesh     = Array{Float32,1}
const FrequencyMeshCplx = Array{Complex{Float32},1}
const SingleParticleFunction = Array{Complex{Float32},3}
const TwoParticleFunction = Array{Complex{Float32},3}  # single orbital for now, depends on 3 frequencies

#new objects ############################################
struct Fockstates
   Nstates    ::Int64
   Nmax       ::Int64  # maximum number of electrons
   norb       ::Int64
   nflavors   ::Int64  # the 'width' of a fockstate 2*norb
   nstatesS ::Array{Int64,1}    # How many different S numbers available for given particle number
   nstatesNS::Array{Int64,2}    # How many states for quantum number N and S are there
   startNS  ::Array{Int64,3}    # starting index for state N,S,i
   SpinNS   ::Array{Int64,2}    # the spin state S for given N,S
   Spin     ::Array{Int64,1}    # the spin state S for Fockstate i
   NelNS    ::Array{Int64,2}    # the number of Electrons for given N,S
   Nel      ::Array{Int64,1}    # the number of Electrons for Fockstate i

   fockstates::Array{FockElement,1} # Store all states in a contiguous array
end

struct Eigenspace
   # Almost identical to Fockstates, is there a better way? (C++ like inheritance?)
   Nstates    ::Int64 
   Nmax       ::Int64  # maximum number of electrons
   norb       ::Int64
   nstatesS   ::Array{Int64,1}    # How many different S numbers available for given particle number
   dimNS      ::Array{Int64,2}    # Dimension of the sub-Eigenspace for quantum number N and S
   dim        ::Array{Int64,1}    # Dimension/Length of i-th Eigenvector
   startNSeval::Array{Int64,3}    # starting index for eval N,S,i
   startNSevec::Array{Int64,3}    # starting index for eigenstate N,S,i
   startevec  ::Array{Int64,1}    # starting index for i-th eigenstate
   SpinNS     ::Array{Int64,2}    # the spin state S for given N,S
   Spin       ::Array{Int64,1}    # the spin state S for Eigenstate i
   NelNS      ::Array{Int64,2}    # the number of electrons for given N,S
   Nel        ::Array{Int64,1}    # the number of electrons for Eigenstate i

   evals::Array{Eigenvalue,1} # Store all states in a contiguous array
   evecs::Array{EigenvectorElem,1} # Store all states in a contiguous array
end

struct Transitions
   Nmax       ::Int64  # maximum number of electrons
   nstatesS   ::Array{Int64,1}    # How many different S numbers available for given particle number
   dimNS      ::Array{Int64,2}    # number of Transitions in the N,S subspace
   startNS    ::Array{Int64,3}    # starting index for the transition N,S,i
   FromTo     ::Array{Int64,2}    # Store all transitions as (from -> to) e.g. FromTo[i,:] = (13,14)
   EvalFromTo ::Array{Eigenvalue,2}    # Store the eigenvalues in the same way as FromTo
   isFromTo   ::Array{Eigenvalue,2}    # Store the is spin numbers in the same way as FromTo
   overlap    ::Array{Complex{Float32},1} # Store all overlaps in a contiguous array
end
# parameter structs #####################################################

struct SimulationParameters
   U      ::Float64        # local Hubbard interaction
   J      ::Float64        # local Hund's coupling
   Up     ::Float64        # interorbital Hubbard interaction=U-2*J
   t      ::Float64        # hopping parameter (positive)
   mu     ::Float64        # chemical potential
   beta   ::Float64        # inverse temperature
   gf_flav::Array{Int64,1} # flavors (orb/spin) for which the Green's function is calculated (full matrix flav x flav)
end

SimulationParameters(;U,J,Up=U-2*J,t,mu,beta,gf_flav) = SimulationParameters(U,J,Up,t,mu,beta,gf_flav)

struct NumericalParameters
   delta            ::Float64   # broadening parameter, we work above the real axis at w+im*delta
   cutoff           ::Float64   # numerical cutoff for things like Boltzmann weights or other things

   NumericalParameters(;delta,cutoff) = new(delta,cutoff)
end

# Outer Constructor for FrequencyMeshes
# init wf with equidistant real frequency grid
# init iwf with matsubaragrid [0,iwmax] and determine number of freq. such that we reach iwmax 

struct FrequencyMeshes
   # Real frequency mesh
   wf::FrequencyMesh   # frequency mesh

   # Fermionic Matsubara frequency mesh
   iwf::FrequencyMesh # frequency mesh

   # Bosonic Matsubara frequency mesh
   ivf::FrequencyMesh # frequency mesh

   FrequencyMeshes(;nw,wmin,wmax,iwmax,beta) = nw>0 ? new( [wmin+n*(wmax-wmin)/nw for n=0:nw-1], 
                                                      [(2*n+1)*pi/beta for n=0:round(Int32, (iwmax*beta/pi-1)/2)-1],
                                                      [(2*n+0)*pi/beta for n=0:round(Int32, (iwmax*beta/pi-1)/2) ] ) :
                                                 throw(ArgumentError("nw=$nw has to be larger than zero!"))
end
