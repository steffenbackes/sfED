"""
    Transition

Fields
-------------
- **`iFromTo`**     :
- **`nFromTo`**     :
- **`sFromTo`**     :
- **`EvalFromTo`**  :
- **`ExpFromTo`**   :
- **`overlap`**     :
"""
struct Transition
   iFromTo   ::Array{Int64,1}    # [index_from, index_to]
   nFromTo   ::Array{Int64,1}    # particle number from-to
   sFromTo   ::Array{Int64,1}    # spin-index from-to
   EvalFromTo::Array{Eigenvalue,1}
   ExpFromTo::Array{Float64,1}
   overlap   ::Complex{Float64}
end

"""
    Transitions

Fields
-------------
- **`transitions`** : n,s,i from
- **`dictFromTo`**  : Dictionary returning the index in the N,S block for specific transiton (from,to)
- **`dictFrom`**    : Dictionary returning the range  of transitions which start from "from" ->any
"""
struct Transitions
   transitions::Array{Array{Array{Transition,1},1},1}
   dictFromTo ::Array{Array{Dict{Tuple{Int64,Int64},Int64},1},1}
   dictFrom   ::Array{Array{Dict{Int64,UnitRange{Int64}},1},1}
end

"""
    Transitions(eigenspace,fockstates,flavors,pNumerics,expCutoff)

Calculate the overlap elements between all Eigenstates acting on c/c^dagger specified by flavors.
The function returns the overlap array and possibleTransitions for cdagger,c for 
the flavor specified (1=orb1,up, 2=orb1,dn, 3=orb2,up, ....)
"""
function Transitions(flavor::Int64,
                     eigenspace::Eigenspace,
                     fockstates::Fockstates,
                     beta::Float64,
                     pNumerics::NumericalParameters,
                     expCutoff::Int64)
	# if flavor>0 we calculate c^dagger
	# if flavor<0 we calculate c

   Nmax = fockstates.norb*2
   nstates = fockstates.Nstates
   E0 = eigenspace.E0
   transitions = [ [] for n=0:Nmax]  # n,s,i
   dictFromTo = [ [] for n=0:Nmax]  # n,s,i
   dictFrom = [ [] for n=0:Nmax]  # n,s,i

   nst = sum([ noSpinConfig(n,Nmax) for n=0:Nmax  ])
   ii = 0
   for n1=0:Nmax
      for s1=1:noSpinConfig(n1,Nmax)
         if ii%(max(1,round(Int64,nst/100.0))) == 0
               print("\r"*lpad(round(Int64,ii*100.0/nst),4)*"%")
         end
         ii +=1
         push!( transitions[n1+1], [] )
         push!( dictFromTo[n1+1], Dict() )
         push!( dictFrom[n1+1], Dict() )
         dS = 2*(abs(flavor)%2)-1 # spin increases by 1 for crea if up flavor, etc
 
         # #  <2|cmat|1> ##############################################
         n2 = n1 + sign(flavor)
         S2 = spinConfig(s1,n1,Nmax) + dS*sign(flavor)
         s2 = indexSpinConfig(S2,n2,Nmax)               # this function returns -1 if state n2,S2 does not exist

         if s2>0
            cmat = getCCdaggerMat(flavor,fockstates.states[n2+1][s2], fockstates.states[n1+1][s1])

            dictFTList = []
            dictFList = []
            noi = 0
            for i=1:length(eigenspace.evals[n1+1][s1])                # now loop over the |1> subspace
               Efrom = eigenspace.evals[n1+1][s1][i] - E0
               expFrom = exp(-beta*Efrom)
               jcounter = 0
               for j=1:length(eigenspace.evals[n2+1][s2])             # now loop over the <2| subspace
                  Eto = eigenspace.evals[n2+1][s2][j] - E0
                  expTo = exp(-beta*Eto)

                  if expCutoff==0 || expTo+expFrom>pNumerics.cutoff
                     ovrlp = dot( eigenspace.evecs[n2+1][s2][j], cmat * eigenspace.evecs[n1+1][s1][i] ) # overlap
                     if abs(ovrlp) > pNumerics.cutoff
                        push!(transitions[n1+1][s1], Transition([i,j],[n1,n2],[s1,s2],[Efrom,Eto],[expFrom,expTo],ovrlp) )     # save transition
                        push!(dictFTList, ( (i,j), length(transitions[n1+1][s1])   ) )
                        jcounter += 1
                        #println(fockstates.states[n1+1][s1][i]," -> ",fockstates.states[n2+1][s2][j], ": ",round.([expFrom,expTo],digits=3)," : ",ovrlp)
                     end # overlap cutoff
                  end # exp cutoff if needed
               end # j states <2|
               if jcounter>0
                  push!(dictFList, ( i, (length(transitions[n1+1][s1])-jcounter+1:length(transitions[n1+1][s1])) ) )
                  jcounter = 0
               end
            end # i states |1>

            dictFromTo[n1+1][s1] = Dict(dictFTList)
            dictFrom[n1+1][s1] = Dict(dictFList)
            #println("N1=",n1,", s1=",s1, ", N2=",n2,", s2=",s2," : ", dictFTList, "  ovrlp: ",transitions[n1+1][s1][1].overlap)
            #println("N1=",n1,", S1=",s1, ", N2=",n2,", S2=",s2," : ", dictFromTo[n1+1][s1])
            #println("N1=",n1,", S1=",s1, ", N2=",n2,", S2=",s2," : ", dictFromt[n1+1][s1])
            #println(" ")
         end # n2states >0
      end # s1
   end # n1
   return Transitions(transitions,dictFromTo,dictFrom)
end


################################################################################
#                               Auxilliary Function                            #
################################################################################

function get1pGFTransitions(flavors::Array{Int64,1},
                     eigenspace::Eigenspace,
                     fockstates::Fockstates,
                     beta::Float64,
                     pNumerics::NumericalParameters)
   # we return c^dagger and c for all flavors
   transitions = Transitions[]
   for m in flavors
      push!( transitions, Transitions( m,eigenspace,fockstates,beta,pNumerics,1) ) #c^dagger
      push!( transitions, Transitions(-m,eigenspace,fockstates,beta,pNumerics,1) ) # c
   end
   println("\rdone!")
   return transitions
end

function get2pGFTransitions(flavor::Int64,
                            eigenspace::Eigenspace,
                            fockstates::Fockstates,
                            beta::Float64,
                            pNumerics::NumericalParameters)

   # we return c^dagger_up, c_up, d^dagger_dn and c_dn for the given flavor, no Boltzmann cutoff
   transitions = Transitions[]
   push!( transitions, Transitions( flavor  ,eigenspace,fockstates,beta,pNumerics,0) ) #c^dagger up
   push!( transitions, Transitions(-flavor  ,eigenspace,fockstates,beta,pNumerics,0) ) # c up
   push!( transitions, Transitions( flavor+1,eigenspace,fockstates,beta,pNumerics,0) ) #c^dagger dn
   push!( transitions, Transitions(-flavor-1,eigenspace,fockstates,beta,pNumerics,0) ) # c dn
   println("\rdone!")
   return transitions
end
