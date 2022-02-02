#typedefs
const FrequencyMesh     = Array{Float64,1}

const FrequencyMeshCplx = Array{Complex{Float64},1}

# parameter structs #####################################################

"""
    ModelParameters

Fields
-------------
- **`U`**       : local Hubbard interaction
- **`J`**       : local Hund's coupling
- **`Up`**      : interorbital Hubbard interaction=U-2*J
- **`t`**       : hopping parameter (positive)
- **`mu`**      : chemical potential
- **`beta`**    : inverse temperature
- **`aim`**     : Anderson IMpurity model: if equal one, apply chemical potential only to first orbital (rest is considered as bath)
- **`gf_flav`** : flavors (orb/spin) for which the Green's function is calculated (full matrix flav x flav)
"""
struct ModelParameters
   U      ::Float64
   J      ::Float64
   Up     ::Float64
   t      ::Float64
   mu     ::Float64
   beta   ::Float64
   aim    ::Int64
   gf_flav::Array{Int64,1}
end

ModelParameters(;U,J,Up=U-2*J,t,mu,beta,aim,gf_flav) = ModelParameters(U,J,Up,t,mu,beta,aim,gf_flav)


"""
    NumericalParameters

Fields
-------------
- **`delta`**  : broadening parameter, we work above the real axis at w+im*delta
- **`cutoff`** : numerical cutoff for things like Boltzmann weights or other things
"""
struct NumericalParameters
   delta            ::Float64
   cutoff           ::Float64

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
                                                      [(2*n+1)*pi/beta for n=0:round(Int64, (iwmax*beta/pi-1)/2)-1],
                                                      [(2*n+0)*pi/beta for n=0:round(Int64, (iwmax*beta/pi-1)/2) ] ) :
                                                 throw(ArgumentError("nw=$nw has to be larger than zero!"))
end
