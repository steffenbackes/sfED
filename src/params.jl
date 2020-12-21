#typedefs
const FrequencyMesh     = Array{Float32,1}
const FrequencyMeshCplx = Array{Complex{Float32},1}

# parameter structs #####################################################

struct SimulationParameters
   U      ::Float64        # local Hubbard interaction
   J      ::Float64        # local Hund's coupling
   Up     ::Float64        # interorbital Hubbard interaction=U-2*J
   t      ::Float64        # hopping parameter (positive)
   mu     ::Float64        # chemical potential
   beta   ::Float64        # inverse temperature
   aim    ::Int64        # Anderson IMpurity model: if equal one, apply chemical potential only to first orbital (rest is considered as bath)
   gf_flav::Array{Int64,1} # flavors (orb/spin) for which the Green's function is calculated (full matrix flav x flav)
end

SimulationParameters(;U,J,Up=U-2*J,t,mu,beta,aim,gf_flav) = SimulationParameters(U,J,Up,t,mu,beta,aim,gf_flav)

struct NumericalParameters
   delta            ::Float32   # broadening parameter, we work above the real axis at w+im*delta
   cutoff           ::Float32   # numerical cutoff for things like Boltzmann weights or other things

   NumericalParameters(;delta,cutoff) = new(delta,cutoff)
end

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
