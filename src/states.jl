
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
   if n>=0 && n<=Nmax
      i = (S + (_nnmax(n,Nmax)))/2 + 1
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


"""
    getNmaxFromAllstates(allstates)

return the maximum number of electrons possible determined from the allstates array
"""
@inline getNmaxFromAllstates(allstates::NSstates) = Int64(size(allstates)[1]-1)

########################################################################

"""
   Fockstates(norb)
Constructor for the Fockstates struct
"""
function Fockstates(; norb )
   norb=norb
   Nstates=4^norb
   Nmax=2*norb
   nflavors=2*norb

   allstates = generateStates(Nmax,Nstates)

   nstatesS = zeros(Int64,Nmax+1)
   nSmax = maximum( length.(allstates) )
   nstatesNS = zeros(Int64,Nmax+1,nSmax)
   SpinNS = zeros(Int64,Nmax+1,nSmax)
   Spin   = zeros(Int64,Nstates)
   NelNS  = zeros(Int64,Nmax+1,nSmax)
   Nel    = zeros(Int64,Nstates)

   dimMax = 0
   ii = 1
   for n=0:Nmax
      nstatesS[n+1] = length(allstates[n+1])
      for is=1:length(allstates[n+1])
         SpinNS[n+1,is] = getSpin(allstates[n+1][is][1])
         NelNS[n+1,is] = sum(allstates[n+1][is][1])
         nstatesNS[n+1,is] = length(allstates[n+1][is])

         if (nstatesNS[n+1,is]>dimMax)
            dimMax = nstatesNS[n+1,is]
         end

         for i=1:nstatesNS[n+1,is]
            Spin[ii] = SpinNS[n+1,is]
            Nel[ii] = NelNS[n+1,is]
            ii += 1
         end

      end
   end

   startNS = zeros(Int64,Nmax+1,nSmax,dimMax)
   fockstates = zeros(FockElement,Nstates*nflavors)

   ii = 1
   for n=0:Nmax
      for is=1:nstatesS[n+1]
         for i=1:nstatesNS[n+1,is]
            istart = ii
            iend = ii+nflavors-1
            fockstates[istart:iend] = copy(allstates[n+1][is][i])

            startNS[n+1,is,i]=istart

            ii+=nflavors
         end
      end
   end
   return Fockstates(Nstates,Nmax,norb,nflavors,nstatesS,nstatesNS,startNS,SpinNS,Spin,NelNS,Nel,fockstates)
end

"""
   getindex(fs,i)
Just return the i-th Fock state in fs
"""
function Base.getindex(fs::Fockstates,i::Int64)
   @assert i>0
   @assert i<=fs.Nstates
   istart=fs.nflavors*(i-1)+1
   iend = istart+fs.nflavors-1
   return fs.fockstates[istart:iend]
end

"""
   getindex(fs,n,is,i)
Return the i-th Fock state in the N=n,S=S(is) subspace
"""
function Base.getindex(fs::Fockstates,n::Int64,is::Int64,i::Int64)
   @assert n>=0  # n is particle number and can be zero
   @assert n<=fs.Nmax
   @assert is>0
   @assert is<=fs.nstatesS[n+1]
   @assert i>0
   @assert i<=fs.nstatesNS[n+1,is]
   istart = fs.startNS[n+1,is,i]
   iend = istart+fs.nflavors-1
   return fs.fockstates[istart:iend]
end

"""
   getindex(fs,n,is,:)
Return the Fock states in the N=n,S=S(is) subspace
"""
function Base.getindex(fs::Fockstates,n::Int64,is::Int64,::Colon)
   @assert n>=0  # n is particle number and can be zero
   @assert n<=fs.Nmax
   @assert is>0
   @assert is<=fs.nstatesS[n+1]
   istart = fs.startNS[n+1,is,1]
   iend = istart+fs.nflavors-1
   return [ fs.fockstates[istart+(i-1)*fs.nflavors:iend+(i-1)*fs.nflavors] for i=1:fs.nstatesNS[n+1,is]]
end

"""
   nstatesNS(fs,N,S)
Return number of states for given N and S
"""
function nstatesNS(fs::Fockstates,N::Int64,S::Int64)
   is = indexSpinConfig(S,N,fs.Nmax)
   if ( is != -1 ) # then it's a valid N,S combination for these states
      return fs.nstatesNS[N+1,is]
   else
      return 0
   end
end


#"""
#   ==(Transition,Transition)
#"""
#Base.:(==)(t1::Transition, t2::Transition) = t1.n1==t2.n1 && t1.n2==t2.n2

"""
   getindex(tr,n,is,i)
Return the i-th transition from N,S -> dN,dS: (indexFrom, indexTo, Efrom,Eto,isfrom, isto, overlap)
"""
function Base.getindex(tr::Transitions,n::Int64,is::Int64,i::Int64)
   @assert n>=0  # n is particle number and can be zero
   @assert n<=tr.Nmax
   @assert is>0
   @assert is<=tr.nstatesS[n+1]
   @assert i>0
   @assert i<=tr.dimNS[n+1,is]
   ii = tr.startNS[n+1,is,i]
   return [tr.FromTo[ii,1], tr.FromTo[ii,2], tr.EvalFromTo[ii,1], tr.EvalFromTo[ii,2], tr.isFromTo[ii,1], tr.isFromTo[ii,2], tr.overlap[ii] ]
end

"""
   getindex(tr,n,is,:)
Return transitions from N,S -> dN,dS
"""
function Base.getindex(tr::Transitions,n::Int64,is::Int64,::Colon)
   @assert n>=0  # n is particle number and can be zero
   @assert n<=tr.Nmax
   @assert is>0
   @assert is<=tr.nstatesS[n+1]
   return [ tr[n,is,i] for i=1:tr.dimNS[n+1,is]]
end

"""
   getindex(tr,n,is,ifrom,ito)
Return the overlap for the transition (ifrom->ito) from N,S -> dN,dS
Return zero if not existing
"""
function Base.getindex(tr::Transitions,n::Int64,is::Int64,ifrom::Int64, ito::Int64)
   @assert n>=0  # n is particle number and can be zero
   @assert n<=tr.Nmax
   @assert is>0
   @assert is<=tr.nstatesS[n+1]
   istart = tr.startNS[n+1,is,1]
   iend = istart+tr.dimNS[n+1,is]-1
   ovrlp = 0
   for it in findall(x->x==[ifrom,ito], [ tr.FromTo[istart:iend,:][t,:] for t=1:iend-istart+1 ])
      ovrlp = tr.overlap[tr.startNS[n+1,is,it]]
   end
   return ovrlp
end

"""
   ntransitionsNS(tr,N,S)
Return number of transitions for given N and S
"""
function ntransitionsNS(tr::Transitions,N::Int64,S::Int64)
   is = indexSpinConfig(S,N,tr.Nmax)
   if ( is != -1 ) # then it's a valid N,S combination for these states
      return tr.dimNS[N+1,is]
   else
      return 0
   end
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
   Nmax=fockstates.Nmax
   norb=fockstates.norb
   nstatesS=copy(fockstates.nstatesS)
   dimNS=copy(fockstates.nstatesNS)
   SpinNS=copy(fockstates.SpinNS)
   Spin=copy(fockstates.Spin)
   NelNS=copy(fockstates.NelNS)
   Nel=copy(fockstates.Nel)

   nSmax = size(fockstates.startNS)[2]
   dimMax = size(fockstates.startNS)[3]
   startNSeval = zeros(Int64,Nmax+1,nSmax,dimMax)
   startNSevec = zeros(Int64,Nmax+1,nSmax,dimMax)
   startevec = zeros(Int64,Nstates)
   dim = zeros(Int64,Nstates)

   #Determine length of 1d eigenvector array
   lenEvecs = 0
   for n=0:Nmax
      for s=1:nstatesS[n+1]
         lenEvecs += dimNS[n+1,s]^2
      end
   end

   # Determine starting points
   ii_eval = 1
   ii_evec = 1
   for n=0:Nmax
      for is=1:nstatesS[n+1]
         for i=1:dimNS[n+1,is]
            startNSeval[n+1,is,i] = ii_eval

            startNSevec[n+1,is,i]=ii_evec
            startevec[ii_eval]=ii_evec

            dim[ii_eval] = dimNS[n+1,is]

            ii_eval += 1
            ii_evec+=dimNS[n+1,is]
         end
      end
   end
   evallist,eveclist = getEvalveclist(eps,tmatrix,Umatrix,Jmatrix,mu,fockstates,lenEvecs,pNumerics)

   return Eigenspace(Nstates,Nmax,norb,nstatesS,dimNS,dim,startNSeval,startNSevec,startevec,SpinNS,Spin,NelNS,Nel,evallist,eveclist)
end

"""
   getindex(es,i)
Just return the i-th eigenvalue and state
"""
function Base.getindex(es::Eigenspace,i::Int64)
   @assert i>0
   @assert i<=es.Nstates
   istart=es.startevec[i]
   iend = istart+es.dim[i]-1
   return es.evals[i], es.evecs[istart:iend]
end

"""
   getindex(es,n,is,i)
Return the i-th Eigenvalue and Eigenvector in the N=n,S=S(is) subspace
"""
function Base.getindex(es::Eigenspace,n::Int64,is::Int64,i::Int64)
   @assert n>=0  # n is particle number and can be zero
   @assert n<=es.Nmax
   @assert is>0
   @assert is<=es.nstatesS[n+1]
   @assert i>0
   @assert i<=es.dimNS[n+1,is]
   istart = es.startNSevec[n+1,is,i]
   iend = istart+es.dimNS[n+1,is]-1

   return es.evals[es.startNSeval[n+1,is,i]], es.evecs[istart:iend]
end

"""
   getindex(es,n,is,:)
Return the Eigenvalues and -states in the N=n,S=S(is) subspace
"""
function Base.getindex(es::Eigenspace,n::Int64,is::Int64,::Colon)
   @assert n>=0  # n is particle number and can be zero
   @assert n<=es.Nmax
   @assert is>0
   @assert is<=es.nstatesS[n+1]
   istart = es.startNSevec[n+1,is,1]
   iend = istart+es.dimNS[n+1,is]-1

   return es.evals[es.startNSeval[n+1,is,1] : es.startNSeval[n+1,is,1]+es.dimNS[n+1,is]-1],
          [ es.evecs[istart+(i-1)*es.dimNS[n+1,is] : iend+(i-1)*es.dimNS[n+1,is]] for i=1:es.dimNS[n+1,is] ]
end

"""
   nstatesNS(es,N,S)
Return number of states for given N and S
"""
function nstatesNS(es::Eigenspace,N::Int64,S::Int64)
   is = indexSpinConfig(S,N,es.Nmax)
   if ( is != -1 ) # then it's a valid N,S combination for these states
      return es.dimNS[N+1,is]
   else
      return 0
   end
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

