struct Transition
   iFromTo   ::Array{Int64,1}    # [index_from, index_to]
   sFromTo   ::Array{Int64,1}    # spin-index from-to
   EvalFromTo::Array{Eigenvalue,1}
   overlap   ::Complex{Float32}
end

struct Transitions
   transitions::Array{Array{Array{Transition,1},1},1}       # n,s,i from
end

## Main functions & constructors ###############################################
"""
    getPossibleTransitions(eigenspace,fockstates,flavors,pNumerics,expCutoff)

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

   for n1=0:Nmax
      #if n1%(max(1,round(Int64,Nmax/100.0))) == 0
      #      print("\r"*lpad(round(Int64,n1*100.0/Nmax),4)*"%")
      #end
      for s1=1:noSpinConfig(n1,Nmax)
         push!( transitions[n1+1], [] )
         dS = 2*(abs(flavor)%2)-1 # spin increases by 1 for crea if up flavor, etc
 
         # #  <2|cmat|1> ##############################################
         n2 = n1 + sign(flavor)
         S2 = spinConfig(s1,n1,Nmax) + dS*sign(flavor)
         s2 = indexSpinConfig(S2,n2,Nmax)               # this function returns -1 if state n2,S2 does not exist

         if s2>0
            cmat = getCCdaggerMat(flavor,fockstates.states[n2+1][s2], fockstates.states[n1+1][s1])

            for i=1:length(eigenspace.evals[n1+1][s1])                # now loop over the |1> subspace
               Efrom = eigenspace.evals[n1+1][s1][i] - E0
               for j=1:length(eigenspace.evals[n2+1][s2])             # now loop over the <2| subspace
                  Eto = eigenspace.evals[n2+1][s2][j] - E0

                  if expCutoff==0 || (exp(-beta*Efrom))+exp(-beta*Eto)>pNumerics.cutoff
                     ovrlp = dot( eigenspace.evecs[n2+1][s2][j], cmat * eigenspace.evecs[n1+1][s1][i] ) # overlap
                     if abs(ovrlp) > pNumerics.cutoff
                        push!(transitions[n1+1][s1], Transition([i,j],[s1,s2],[Efrom,Eto],ovrlp) )     # save transition
                     end # overlap cutoff
                  end # exp cutoff if needed
               end # j states <2|
            end # i states |1>
         end # n2states >0
      end # s1
   end # n1

   return Transitions(transitions)
end

#########################################################################################
function get1pGFTransitions(flavors::Array{Int64,1},
                     eigenspace::Eigenspace,
                     fockstates::Fockstates,
                     beta::Float64,
                     pNumerics::NumericalParameters)
	# if flavor>0 we calculate c^dagger
	# if flavor<0 we calculate c

   # we return c^dagger and c for all flavors
   transitions = Transitions[]
   for m in flavors
      push!( transitions, Transitions( m,eigenspace,fockstates,beta,pNumerics,1) ) #c^dagger
      push!( transitions, Transitions(-m,eigenspace,fockstates,beta,pNumerics,1) ) # c
   end
   return transitions
end


#################
