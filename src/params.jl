#typedefs
const CAmatrix    = SparseMatrixCSC{Int8,Int64}   # matrix type for creation/annihilation matrix
const Hamiltonian = SparseMatrixCSC{Complex{Float32},Int64}   # matrix type for Hamiltonian matrix

const FockElement = Int8                                  # one number element of a Fock state
const Fockstate   = Array{FockElement,1}                  # The standard Fock state e.g. [00101011] 
const NSstates    = Array{Array{Array{Fockstate,1},1},1}  # An array of arrays for all Fock states for Quantum numbers N,S: typeof( allstates[N][S][i] )=Fockstate

const Eigenvalue        = Float32             
const Eigenvector       = Array{Complex{Float32},1}
const EigenvectorMatrix = Array{Complex{Float32},2}

const FrequencyMesh     = Array{Float32,1}
const FrequencyMeshCplx = Array{Complex{Float32},1}
const SingleParticleFunction = Array{Complex{Float32},3}

struct ModelParameters
	norb   ::UInt32    # Number of Orbitals
	Nmax   ::UInt32    # Maximal number of electrons=2*norb
	Nstates::UInt64     # Total number of possible Fock states=4^norb  (0,up,dn,updn for each state)

	#Constructor
	ModelParameters(;norb) = new(norb, 2*norb, 4^norb)
end

struct SimulationParameters
	U      ::Float64        # local Hubbard interaction
	J      ::Float64        # local Hund's coupling
	Up     ::Float64        # interorbital Hubbard interaction=U-2*J
	t      ::Float64        # hopping parameter (positive)
	mu     ::Float64        # chemical potential
	beta   ::Float64        # inverse temperature
	gf_orbs::Array{Int64,1} # orbitals for which the Green's function is calculated (full matrix norb x norb)
	
	# Constructor
	SimulationParameters(;U,J,t,mu,beta,gf_orbs) = new(U,J,U-2*J,t,mu,beta,gf_orbs)
end

struct NumericalParameters
	delta            ::Float64   # broadening parameter, we work above the real axis at w+im*delta
	cutoff           ::Float64   # numerical cutoff for things like Boltzmann weights or other things
	nevalsPerSubspace::UInt64     # how much eigenvalues to obtain for each subspace (affects ARPACK directly)
	nevalsTotalMax   ::UInt64     # how much eigenvalues to keep in total from all subspaces (affects mostly memory and Green's function)

	NumericalParameters(;delta,cutoff,nevalsPerSubspace,nevalsTotalMax) = new(delta,cutoff,nevalsPerSubspace,nevalsTotalMax)
end

# Outer Constructor for FrequencyMeshes
# init wf with equidistant real frequency grid
# init iwf with matsubaragrid [0,iwmax] and determine number of freq. such that we reach iwmax 

struct FrequencyMeshes
	# Real frequency mesh
	wf::FrequencyMesh   # frequency mesh

	# Fermionic Matsubara frequency mesh
	iwf::FrequencyMesh # frequency mesh

	FrequencyMeshes(;nw,wmin,wmax,iwmax,beta) = new( [wmin+n*(wmax-wmin)/nw for n=0:nw-1], 
	                                                [(2*n+1)*pi/beta for n=0:round(Int32, (iwmax*beta/pi-1)/2)-1] )
end
