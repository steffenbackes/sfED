
struct ModelParameters
	norb   ::Int64     # Number of Orbitals
	Nmax   ::Int64     # Maximal number of electrons=2*norb
	Nstates::Int64     # Total number of possible Fock states=4^norb  (0,up,dn,updn for each state)

	#Constructor
	ModelParameters(;norb) = new(norb, 2*norb, 4^norb)
end

struct SimulationParameters
	U     ::Float64        # local Hubbard interaction
	J     ::Float64        # local Hund's coupling
	Up    ::Float64        # interorbital Hubbard interaction=U-2*J
	t     ::Float64        # hopping parameter (positive)
	mu    ::Float64        # chemical potential
	beta  ::Float64        # inverse temperature
	
	# Constructor
	SimulationParameters(;U,J,t,mu,beta) = new(U,J,U-2*J,t,mu,beta)
end

struct NumericalParameters
	delta            ::Float64   # broadening parameter, we work above the real axis at w+im*delta
	cutoff           ::Float64   # numerical cutoff for things like Boltzmann weights or other things
	nevalsPerSubspace::Int64     # how much eigenvalues to obtain for each subspace (affects ARPACK directly)
	nevalsTotalMax   ::Int64     # how much eigenvalues to keep in total from all subspaces (affects mostly memory and Green's function)

	NumericalParameters(;delta,cutoff,nevalsPerSubspace,nevalsTotalMax) = new(delta,cutoff,nevalsPerSubspace,nevalsTotalMax)
end

# Outer Constructor for FrequencyMeshes
# init wf with equidistant real frequency grid
# init iwf with matsubaragrid [0,iwmax] and determine number of freq. such that we reach iwmax 

struct FrequencyMeshes
	# Real frequency mesh
	wf::Array{Float64,1}   # frequency mesh

	# Fermionic Matsubara frequency mesh
	iwf::Array{Float64,1} # frequency mesh

	FrequencyMeshes(;nw,wmin,wmax,iwmax,beta) = new( [wmin+n*(wmax-wmin)/nw for n=0:nw-1], 
	                                                [(2*n+1)*pi/beta for n=0:round(Int64, (iwmax*beta/pi-1)/2)-1] )
end
