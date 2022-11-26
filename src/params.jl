# typedefs ##############################################################
const FrequencyMesh     = Vector{Float64}
const FrequencyMeshCplx = Vector{ComplexF64}

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
end
NumericalParameters(;delta,cutoff) = NumericalParameters(delta,cutoff)


"""
    FrequencyMeshes

Frequency meshes for 1 and 2 particle calculations. Constructed via `FrequencyMeshes(;nw, wmin, wmax, iwmax, beta)` with `nw` being the length of the real frequency mesh starting at `wmin` and ending at `wmax`.
    `iwmax` gives the number of negative frequencies (there are equally many positive frequencies generated in the bosonic case and `iwmax-1` positive frequencies in the fermionic one) for the fermionic and bosonic mesh at temperature `beta`.

Fields
-------------

- **`wf`**  : Real frequency mesh
- **`iwf`** : Fermionic Matsubara frequency mesh
- **`ivf`** : Bosonic Matsubara frequency mesh
"""
struct FrequencyMeshes
    wf::FrequencyMesh
    iwf::FrequencyMesh
    ivf::FrequencyMesh
    function FrequencyMeshes(;nw, wmin, wmax, iwmax, beta) 
        if nw > 0
            new([wmin+n*(wmax-wmin)/nw for n=0:nw-1], 
                [(2*n+1)*pi/beta for n=-iwmax:iwmax],
                [(2*n+0)*pi/beta for n=-iwmax:iwmax-1] )
        else
            throw(ArgumentError("nw=$nw has to be larger than zero!"))
        end
    end
end
