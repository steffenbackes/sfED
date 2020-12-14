
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
indexSpinConfig(S::Int64,n::Int64,Nmax::Int64)::Int64 = round(Int64,  (S + (_nnmax(n,Nmax)))/2 + 1)


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


"""
    getNmaxFromAllstates(allstates)

return the maximum number of electrons possible determined from the allstates array
"""
@inline getNmaxFromAllstates(allstates::NSstates) = Int64(size(allstates)[1]-1)

########################################################################

########################################################################
# Here we create the array that stores all states as an integer array sorted by N,S quantum numbers
function generateStates(pModel::ModelParameters)

   Nmax = pModel.Nmax
   Nstates = pModel.Nstates

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
########################################################################



#######################################################################
# Return permutation that sorts the eigenvalue/vector array according to N,S quantum numbers (increasing order), then energy
function getNSperm(evallist::Array{Array{Eigenvalue,1},1})
   NS=[]
   for i=1:length(evallist)
      push!(NS, [evallist[i][2],evallist[i][3],evallist[i][1]] )   # extract the n,s,E quantum number into an array
   end
   return sortperm(NS)
end
#######################################################################

