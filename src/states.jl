const FockElement = Int8                                  # one number element of a Fock state
const Fockstate   = Array{FockElement,1}                  # The standard Fock state e.g. [00101011] 
const NSstates    = Array{Array{Array{Fockstate,1},1},1}  # An array of arrays for all Fock states for Quantum numbers N,S: typeof( allstates[N][S][i] )=Fockstate

const Eigenvalue        = Float32             
const EigenvectorElem   = Complex{Float32}
const Eigenvector       = Array{EigenvectorElem,1}
const EigenvectorMatrix = Array{EigenvectorElem,2}

struct Fockstates
   Nstates    ::Int64
   norb       ::Int64

   states::NSstates   # states[n][s][i]  (electron number, spin index, i-th state in NS block)
end

struct Eigenspace
   Nstates    ::Int64 
   norb       ::Int64
   E0         ::Eigenvalue

   evals::Array{Array{Array{Eigenvalue,1},1},1}       # n,s,i
   evecs::Array{Array{Array{Eigenvector,1},1},1}      # n,s,i
end

# ============================ Comments ============================
#   
#
#
# ======================= Auxilliary Function =======================

function _nnmax(n::Int64, Nmax::Int64)::Int64
   n >= 0 || throw(DomainError(n, "argument n must be nonnegative")) 
   Nmax >= 0 || throw(DomainError(Nmax, "argument Nmax must be nonnegative")) 
   n <= Nmax || throw(DomainError(n-Nmax, "n must be smaller than Nmax")) 
   return Nmax/2-abs(n-Nmax/2)
end

# ========================= Main Functions =========================

"""
    noSpinConfig(n, Nmax)

Number of configurations with different total spin S possible for given `n` and `Nmax`.

# Examples
```
julia> sfED.noSpinConfig(2,4)
3
```
"""
function noSpinConfig(n::Int64, Nmax::Int64)::Int64
    return round(Int64, _nnmax(n,Nmax)+1 )
end

"""
    spinConfig(i,n,Nmax)

`i`-th possible total spin value S for given `n`, `nMax`. `i` starts at 1

# Examples
```
julia> sfED.spinConfig(1,2,4)
-2
julia> sfED.spinConfig(2,2,4)
0
julia> sfED.spinConfig(3,2,4)
2
```
"""
function spinConfig(i::Int64,n::Int64,Nmax::Int64)::Int64
   nTot = noSpinConfig(n,Nmax)
   (i > 0 && i <= nTot) || throw(DomainError(nTot-i, "i must be larger than zero and "*
                                             "smaller than the number of spin configurations.")) 
   return round(Int64, -(_nnmax(n,Nmax)) + (i-1)*2 )
end

"""
    indexSpinConfig(S,n,Nmax)

Index `i` corresponding to a spin config `S` for given `n` and `nMax`

# Example
```
julia> indexSpinConfig(spinConfig(2,2,4),2,4)
2
```
"""
function indexSpinConfig(S::Int64,n::Int64,Nmax::Int64)::Int64
   i = -1
   if n>=0 && n<=Nmax # we don't want an assert, if n is bad return -1
      i = (S + (Nmax/2-abs(n-Nmax/2)))/2 + 1
      # this if is very ugly, any recommendation?
      if abs(round(Int64,i) - i  ) > 0.001  || i<1 || i>noSpinConfig(n,Nmax)  # if one requests a spin state that does not exit return -1, we need this behavior when summing over the N-1,S-1 states
         i=-1
      end
   end
   return i
end


"""
    getSpin(states)

total spin `S` of `state`.
"""
@inline getSpin(state::Fockstate) = sum((2*(i%2)-1)*state[i] for i in 1:length(state))


"""
    getCsign(i,state)

return the sign when creating/annihilating an electron at position `i` in `state`.
"""
@inline getCsign(i::Int64,state::Fockstate) = 1-2*(sum(state[1:i-1])%2)

#####################################################################

"""
   Fockstates(norb)
Constructor for the Fockstates struct
"""
function Fockstates(; norb )
   norb=norb
   Nstates=4^norb

   Nmax=2*norb
   allstates = generateStates(Nmax,Nstates)

   return Fockstates(Nstates,norb,allstates)
end

#####################################################################

"""
   Eigenspace(eps,tmatrix,Umatrix,Jmatrix,pSimulation.mu,fockstates,pNumerics)
Constructor for the Eigenspace struct
"""
function Eigenspace(eps::Array{Float64,1},tmatrix::Array{Float64,2},
                    Umatrix::Array{Float64,2},Jmatrix::Array{Float64,2},
                    mu::Float64,
                    fockstates::Fockstates,
                    pNumerics::NumericalParameters)
   Nstates=fockstates.Nstates
   norb=fockstates.norb

   evals,evecs = getEvalveclist(eps,tmatrix,Umatrix,Jmatrix,mu,fockstates,pNumerics)
   E0 = minimum([ e for nse in evals for se in nse for e in se  ])

   return Eigenspace(Nstates,norb,E0,evals,evecs)
end



########################################################################
# Here we create the array that stores all states as an integer array sorted by N,S quantum numbers
function generateStates(Nmax::Int64,Nstates::Int64)
   # Create an empty list that will store all states, sorted for given particle number and spin
   allstates::NSstates = []
   for n=0:Nmax
      push!( allstates, [] )
      for j=1:noSpinConfig(n,Nmax)
         push!( allstates[n+1], [] )
      end
   end # n
   
   # now create integer arrays for the states and store them according to N,S 
   #(i,e, [0,1,0,1] has N=2 and S=-2)
   # the 1,3,5,7, ... elements are up electrons
   # the 2,4,6,8, ... elements are dn electrons
   # 1,2 is first orbital, 3,4 second, etc....
   for i=1:Nstates
      state = parse.(FockElement, split( bitstring(i-1)[end-Nmax+1:end],"") )
      N = sum(state)
      S = getSpin(state)
      iS = indexSpinConfig(S,N,Nmax)
      push!( allstates[N+1][iS], state )
   end
   return allstates
end



