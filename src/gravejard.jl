
"""
    getGF(evallist, eveclist, allstates, pModel, pSimulation, pFreq, pNumerics)

Evaluation of the Lehmann representation for the interacting finite-temperature Green's function
We sum over all Eigenstates `evallist` with electron number N and all other eigenstates n2 with electron number N-1
Since we have spin degeneracy, we only calculate the up GF, so we annihilate one up spin, so we also only sum over all S-1 states for given S(n1)
"""
function getGFnonoptim( evallist::Array{Array{Eigenvalue,1},1},
                eveclist::Array{Eigenvector,1},
               allstates::NSstates,
               pModel::ModelParameters,
               pSimulation::SimulationParameters,
               pFreq::FrequencyMeshes,
               pNumerics::NumericalParameters)::Tuple{SingleParticleFunction,SingleParticleFunction,Array{Float64,2} }

   norb=length(pSimulation.gf_flav)
   Nmax=pModel.Nmax
   beta = pSimulation.beta

   NSvalues   = Int64[-1,-1,-1,-1]
   cmatrix    = CAmatrix[]
   cdagmatrix = CAmatrix[]
   gfdiagnorm                    = zeros(Float64,norb)                   # check normalization of GF
   gf_w::SingleParticleFunction  = zeros(norb,norb,length(pFreq.wf))
   gf_iw::SingleParticleFunction = zeros(norb,norb,length(pFreq.iwf))

   #evalContributions::Array{Array{Float64,1},1} = []
   evalContributions = Array{Float64}(undef, (length(evallist), 4) )

   # First create a sortting for all Eigenstates in order of increasing N,S, which makes the GF generation more efficient
   NSperm = getNSperm(evallist)

   E0 = minimum( first.(evallist) )

   for n1=1:length(evallist)
       if n1%(round(Int64,length(evallist)/10)) == 0
            print("\r"*lpad(round(Int64,n1*100.0/length(evallist)),4)*"%")
      end

      E1    = evallist[NSperm[n1]][1] -E0     # shift by E0 to avoid overflow
      N1    = round(Int64,evallist[NSperm[n1]][2])
      s1    = round(Int64,evallist[NSperm[n1]][3])
      evec1 = eveclist[NSperm[n1]]
      S1    = spinConfig(s1,N1,Nmax)
      evalContributions[n1,:] = [N1,S1,E1,0.0]
   
      for n2=1:length(evallist)
         E2    = evallist[NSperm[n2]][1] -E0
         N2    = round(Int64,evallist[NSperm[n2]][2])   
         s2    = round(Int64,evallist[NSperm[n2]][3])
         evec2 = eveclist[NSperm[n2]]
         S2    = spinConfig(s2,N2,Nmax)

         expFac = exp(-beta*E1)+exp(-beta*E2)
   
         # Exclude transitions too high in energy or all Eigenstates 2 which are not N-1,S-1 
         if ( expFac > pNumerics.cutoff && N2==N1-1 && S2==S1-1 )
   
            # If we have not been dealing with this N,S combination in the loop before, we need to generate the right Cmatrix and Cdagmatrix
            if ( NSvalues!=[N1,N2,S1,S2] )
               NSvalues = [N1,N2,S1,S2]
   
               # we need to generate cmatrix and cdagmatrix for all orbitals
               dim1 = length(allstates[N1+1][s1])
               dim2 = length(allstates[N2+1][s2])
               cmatrix    = CAmatrix[]   # go from subspace 1 to 2 by annihilation
               cdagmatrix = CAmatrix[]   # go from subspace 2 to 1 by creation
               for m1=1:norb
                  c1 = pSimulation.gf_flav[m1]
                  a1 = pSimulation.gf_flav[m1]
                  push!( cmatrix,       getCmatrix(a1, allstates[N2+1][s2], allstates[N1+1][s1]) )
                  push!( cdagmatrix, getCdagmatrix(c1, allstates[N1+1][s1], allstates[N2+1][s2]) )
               end # m1 loop
            end # if we need to update cmatrix,cdagmatrix
   
            #Now we have the proper N,N-1 state and the c and c^dagger matrices, so evaluate the Lehmann representation
            # Obtain all matrix elements
            for m1=1:norb
               # first the diagonal elements
               evec1cdagevec2 = dot( evec1, cdagmatrix[m1]*evec2 )  # we can reuse this result

               # if this overlap is too small, it doesn't matter what the other is since they get multiplied
               if abs(evec1cdagevec2)>pNumerics.cutoff
   
                  overlap =  evec1cdagevec2 * dot( evec2, cmatrix[m1]*evec1 ) * expFac # include the Boltzmann terms
                  if abs(overlap)>pNumerics.cutoff
                     gf_w[m1,m1,:]  += overlap ./ ( pFreq.wf     .+ (pNumerics.delta*im - E1 + E2) )
                     gf_iw[m1,m1,:] += overlap ./ ( im*pFreq.iwf .+ (                   - E1 + E2) )
                     evalContributions[n1,4] += abs(overlap)
                     evalContributions[n2,4] += abs(overlap)
                     gfdiagnorm[m1] += real(overlap)
                  end
      
                  # then upper offdiagonal
                  for m2=m1+1:norb
                     overlap =  evec1cdagevec2 * dot( evec2, cmatrix[m2]*evec1 ) * expFac # include the Boltzmann terms
                     if abs(overlap)>pNumerics.cutoff
                        gf_w[m1,m2,:]  += overlap ./ ( pFreq.wf     .+ (pNumerics.delta*im - E1 + E2) )
                        gf_iw[m1,m2,:] += overlap ./ ( im*pFreq.iwf .+ (                   - E1 + E2) )
                        evalContributions[n1,4] += abs(overlap)
                        evalContributions[n1,4] += abs(overlap)
                     end
                  end # m2 loop

               end # if evec1cdagevec2>pNumerics.cutoff
            end # m1 loop
   
         end # if exp(-beta*E) > cutoff
      end # n2 loop
   end # n1 loop
    println("\rdone!")
   
   # copy the offdiagonals (not true copy in julia, just referencing but this is ok)
   for m1=1:norb
      for m2=m1+1:norb
         gf_w[m2,m1,:] = gf_w[m1,m2,:]
         gf_iw[m2,m1,:] = gf_iw[m1,m2,:]
      end
   end
   
   #normalize
   gfdiagnorm ./= getZ(evallist,beta)
   gf_w ./= getZ(evallist,beta)
   gf_iw ./= getZ(evallist,beta)

   # normalize the Matsubara GF tail so that the Dyson equation works for the Selfenergy when the normalization is slightly too small
   for m1=1:norb
      gf_iw[m1,m1,:] ./= gfdiagnorm[m1]
   end

   writeGFnorm(gfdiagnorm, maximum( first.(evallist) )-E0)

   return gf_w, gf_iw, evalContributions
end 

####################################################################################################

"""
    getGF(evallist, eveclist, allstates, pModel, pSimulation, pFreq, pNumerics)

Evaluation of the Lehmann representation for the interacting finite-temperature Green's function
We sum over all Eigenstates `evallist` with electron number N and all other eigenstates n2 with electron number N-1
Since we have spin degeneracy, we only calculate the up GF, so we annihilate one up spin, so we also only sum over all S-1 states for given S(n1)
"""
function getGFNSoptim( evallist::Array{Array{Eigenvalue,1},1},
                eveclist::Array{Eigenvector,1},
               fockstates::Fockstates,
               pSimulation::SimulationParameters,
               pFreq::FrequencyMeshes,
               pNumerics::NumericalParameters)::Tuple{SingleParticleFunction,SingleParticleFunction,Array{Float64,2} }

   norb=length(pSimulation.gf_flav)
   Nmax=fockstates.Nmax
   beta = pSimulation.beta

   NSvalues   = Int64[-1,-1,-1,-1]
   cmatrix    = CAmatrix[]
   cdagmatrix = CAmatrix[]
   gfdiagnorm                    = zeros(Float64,norb)                   # check normalization of GF
   gf_w::SingleParticleFunction  = zeros(norb,norb,length(pFreq.wf))
   gf_iw::SingleParticleFunction = zeros(norb,norb,length(pFreq.iwf))

   #evalContributions::Array{Array{Float64,1},1} = []
   evalContributions = Array{Float64}(undef, (length(evallist), 4) )

   # First create a sortting for all Eigenstates in order of increasing N,S, which makes the GF generation more efficient
   NSperm = getNSperm(evallist)

   E0 = minimum( first.(evallist) )

   #####################################
   NSelements = zeros(Int64,Nmax+1,noSpinConfig(Int64(Nmax/2),Nmax))
   for n1=1:length(evallist)
      E1    = evallist[NSperm[n1]][1] -E0     # shift by E0 to avoid overflow
      N1    = round(Int64,evallist[NSperm[n1]][2])
      s1    = round(Int64,evallist[NSperm[n1]][3])
      S1    = spinConfig(s1,N1,Nmax)

      NSelements[N1+1, s1] += 1
   end
   NSstartpoints = zeros(Int64,Nmax+1,noSpinConfig(Int64(Nmax/2),Nmax))
   NSstartpoints[1,1] = 1
   start=2
   for n=1:Nmax
      for s=1:noSpinConfig(n,Nmax)
         NSstartpoints[n+1,s] = start
         start += NSelements[n+1,s]
      end
   end

   #####################################

   for n1=2:length(evallist)
       if n1%(round(Int64,length(evallist)/10)) == 0
            print("\r"*lpad(round(Int64,n1*100.0/length(evallist)),4)*"%")
      end

      E1    = evallist[NSperm[n1]][1] -E0     # shift by E0 to avoid overflow
      N1    = round(Int64,evallist[NSperm[n1]][2])
      s1    = round(Int64,evallist[NSperm[n1]][3])
      evec1 = eveclist[NSperm[n1]]
      S1    = spinConfig(s1,N1,Nmax)
      evalContributions[n1,:] = [N1,S1,E1,0.0]
   
      n2elem = 0
      if ( indexSpinConfig(S1-1,N1-1,Nmax)>0 )
         n2elem = NSelements[N1+1-1,indexSpinConfig(S1-1,N1-1,Nmax)]
      end

      for in2=1:n2elem
         n2 = NSstartpoints[N1+1-1,indexSpinConfig(S1-1,N1-1,Nmax)]-1+in2
         E2    = evallist[NSperm[n2]][1] -E0
         N2    = round(Int64,evallist[NSperm[n2]][2])   
         s2    = round(Int64,evallist[NSperm[n2]][3])
         evec2 = eveclist[NSperm[n2]]
         S2    = spinConfig(s2,N2,Nmax)

         expFac = exp(-beta*E1)+exp(-beta*E2)
   
         # Exclude transitions too high in energy or all Eigenstates 2 which are not N-1,S-1 
         if ( expFac > pNumerics.cutoff )
   
            # If we have not been dealing with this N,S combination in the loop before, we need to generate the right Cmatrix and Cdagmatrix
            if ( NSvalues!=[N1,N2,S1,S2] )
               NSvalues = [N1,N2,S1,S2]
   
               # we need to generate cmatrix and cdagmatrix for all orbitals
               dim1 = fockstates.nstatesNS[N1+1,s1]
               dim2 = fockstates.nstatesNS[N1+1,s2]
               cmatrix    = CAmatrix[]   # go from subspace 1 to 2 by annihilation
               cdagmatrix = CAmatrix[]   # go from subspace 2 to 1 by creation
               for m1=1:norb
                  c1 = pSimulation.gf_flav[m1]
                  a1 = pSimulation.gf_flav[m1]
                  push!( cmatrix,       getCmatrix(a1, fockstates, N2,s2,N1,s1) )
                  push!( cdagmatrix, getCdagmatrix(c1, fockstates, N1,s1,N2,s2) ) 
               end # m1 loop
            end # if we need to update cmatrix,cdagmatrix
   
            #Now we have the proper N,N-1 state and the c and c^dagger matrices, so evaluate the Lehmann representation
            # Obtain all matrix elements
            for m1=1:norb
               # first the diagonal elements
               evec1cdagevec2 = dot( evec1, cdagmatrix[m1]*evec2 )  # we can reuse this result

               # if this overlap is too small, it doesn't matter what the other is since they get multiplied
               if abs(evec1cdagevec2)>pNumerics.cutoff
   
                  overlap =  evec1cdagevec2 * dot( evec2, cmatrix[m1]*evec1 ) * expFac # include the Boltzmann terms
                  if abs(overlap)>pNumerics.cutoff
                     gf_w[m1,m1,:]  += overlap ./ ( pFreq.wf     .+ (pNumerics.delta*im - E1 + E2) )
                     gf_iw[m1,m1,:] += overlap ./ ( im*pFreq.iwf .+ (                   - E1 + E2) )
                     evalContributions[n1,4] += abs(overlap)
                     evalContributions[n2,4] += abs(overlap)
                     gfdiagnorm[m1] += real(overlap)
                  end
      
                  # then upper offdiagonal
                  for m2=m1+1:norb
                     overlap =  evec1cdagevec2 * dot( evec2, cmatrix[m2]*evec1 ) * expFac # include the Boltzmann terms
                     if abs(overlap)>pNumerics.cutoff
                        gf_w[m1,m2,:]  += overlap ./ ( pFreq.wf     .+ (pNumerics.delta*im - E1 + E2) )
                        gf_iw[m1,m2,:] += overlap ./ ( im*pFreq.iwf .+ (                   - E1 + E2) )
                        evalContributions[n1,4] += abs(overlap)
                        evalContributions[n1,4] += abs(overlap)
                     end
                  end # m2 loop

               end # if evec1cdagevec2>pNumerics.cutoff
            end # m1 loop
   
         end # if exp(-beta*E) > cutoff
      end # n2 loop
   end # n1 loop
    println("\rdone!")
   
   # copy the offdiagonals (not true copy in julia, just referencing but this is ok)
   for m1=1:norb
      for m2=m1+1:norb
         gf_w[m2,m1,:] = gf_w[m1,m2,:]
         gf_iw[m2,m1,:] = gf_iw[m1,m2,:]
      end
   end
   
   #normalize
   gfdiagnorm ./= getZ(evallist,beta)
   gf_w ./= getZ(evallist,beta)
   gf_iw ./= getZ(evallist,beta)

   # normalize the Matsubara GF tail so that the Dyson equation works for the Selfenergy when the normalization is slightly too small
   for m1=1:norb
      gf_iw[m1,m1,:] ./= gfdiagnorm[m1]
   end

   writeGFnorm(gfdiagnorm, maximum( first.(evallist) )-E0)

   return gf_w, gf_iw, evalContributions
end 

####################################################################################################
"""
    getGF(evallist, eveclist, allstates, pModel, pSimulation, pFreq, pNumerics)

Evaluation of the Lehmann representation for the interacting finite-temperature Green's function
We sum over all Eigenstates `evallist` with electron number N and all other eigenstates n2 with electron number N-1
Since we have spin degeneracy, we only calculate the up GF, so we annihilate one up spin, so we also only sum over all S-1 states for given S(n1)
"""
function getGFhalfoptim( evallist::Array{Array{Eigenvalue,1},1},
               overlaps::Array{Complex{Float32},3},
               possibleTransitions::Array{Array{Array{Int64,1},1},1}, 
               NSperm::Array{Int64,1},
               pModel::ModelParameters,
               pSimulation::SimulationParameters,
               pFreq::FrequencyMeshes,
               pNumerics::NumericalParameters)::Tuple{SingleParticleFunction,SingleParticleFunction,Array{Float64,2}}

   nflav=length(pSimulation.gf_flav)
   Nmax=pModel.Nmax
   beta = pSimulation.beta

   gfdiagnorm                    = zeros(Float64,nflav)                   # check normalization of GF
   gf_w::SingleParticleFunction  = zeros(nflav,nflav,length(pFreq.wf))
   gf_iw::SingleParticleFunction = zeros(nflav,nflav,length(pFreq.iwf))

   #evalContributions::Array{Array{Float64,1},1} = []

   E0 = minimum( first.(evallist) )
   evalContributions = Array{Float64}(undef, (length(evallist), 4) )

   # prefill evalContribution array
   for n1=1:length(evallist)
      E1    = evallist[NSperm[n1]][1] -E0 
      N1    = round(Int64,evallist[NSperm[n1]][2])
      s1    = round(Int64,evallist[NSperm[n1]][3])
      S1    = spinConfig(s1,N1,Nmax)
      #push!( evalContributions, [N1,S1,E1,0.0] )
      evalContributions[n1, :] .= [N1,S1,E1,0.0]
   end

   for n1=1:length(evallist)
       if n1%(round(Int64,length(evallist)/10)) == 0
            print("\r"*lpad(round(Int64,n1*100.0/length(evallist)),4)*"%")
      end
      E1    = evallist[NSperm[n1]][1] -E0     # shift by E0 to avoid overflow
      N1    = round(Int64,evallist[NSperm[n1]][2])
      s1    = round(Int64,evallist[NSperm[n1]][3])
      S1    = spinConfig(s1,N1,Nmax)
   
      # now we loop over the flavors
      # GF is defined as in the pdf docu:
      # <n1| c_a |n2><n2| cdag_b |n1>
      for a=1:length(pSimulation.gf_flav)
         # possibleTransitions[2*b-1][n1] contains all states that have overlap with cdag_b |n1>
         # possibleTransitions[2*a-1][n1] contains all states that have overlap with cdag_a |n1>, which is the same as <n1| c_a
         # so we only need to sum over states which are in the intersection of these two sets in principle
         # It actually turns out that performing the restricted n2 sum first and then checking for the b overlap is the fastest,
         # probably because the expFac cutoff helps a lot

         for n2 in possibleTransitions[2*a-1][n1]
            E2    = evallist[NSperm[n2]][1] -E0
            N2    = round(Int64,evallist[NSperm[n2]][2])   
            s2    = round(Int64,evallist[NSperm[n2]][3])
            S2    = spinConfig(s2,N2,Nmax)

            expFac = exp(-beta*E1)+exp(-beta*E2)
   
            if ( expFac > pNumerics.cutoff )
               for b=a:length(pSimulation.gf_flav) # only upper triangular part
                  if (abs(expFac*overlaps[2*b-1,n2,n1])>pNumerics.cutoff)

                     ovrlp = expFac*overlaps[2*a-0,n1,n2]*overlaps[2*b-1,n2,n1]            # <n1|c_a|n2> * <n2|cdag_b|n1>
                     gf_w[a,b,:]  += ovrlp ./ ( pFreq.wf     .+ (pNumerics.delta*im + E1 - E2) )
                     gf_iw[a,b,:] += ovrlp ./ ( im*pFreq.iwf .+ (                   + E1 - E2) )
                     evalContributions[n1,4] += abs(ovrlp)
                     evalContributions[n2,4] += abs(ovrlp)
                     if a==b
                        gfdiagnorm[a] += real(ovrlp)
                     end

                  end # cdag_b transition cutoff
               end # b loop
            end # expFac cutoff

         end # n2 loop
      end # a loop
   end # n1 loop
    println("\rdone!")
   
   # copy the offdiagonals (not true copy in julia, just referencing but this is ok)
   for m1=1:nflav
      for m2=m1+1:nflav
         gf_w[m2,m1,:] = gf_w[m1,m2,:]
         gf_iw[m2,m1,:] = gf_iw[m1,m2,:]
      end
   end
   
   #normalize
   gfdiagnorm ./= getZ(evallist,beta)
   gf_w ./= getZ(evallist,beta)
   gf_iw ./= getZ(evallist,beta)

   # normalize the Matsubara GF tail so that the Dyson equation works for the Selfenergy when the normalization is slightly too small
#   for m1=1:nflav
#      gf_iw[m1,m1,:] ./= gfdiagnorm[m1]
#   end

   writeGFnorm(gfdiagnorm, maximum( first.(evallist) )-E0)

   return gf_w, gf_iw, evalContributions
end 
#########################################################################################
"""
    getGF(evallist, eveclist, allstates, pModel, pSimulation, pFreq, pNumerics)

Evaluation of the Lehmann representation for the interacting finite-temperature Green's function
We sum over all Eigenstates `evallist` with electron number N and all other eigenstates n2 with electron number N-1
Since we have spin degeneracy, we only calculate the up GF, so we annihilate one up spin, so we also only sum over all S-1 states for given S(n1)
"""
function getGFold( evallist::Array{Array{Eigenvalue,1},1},
               overlaps::Array{Complex{Float32},3},
               possibleTransitions::Array{Array{Array{Int64,1},1},1}, 
               NSperm::Array{Int64,1},
               pModel::ModelParameters,
               pSimulation::SimulationParameters,
               pFreq::FrequencyMeshes,
               pNumerics::NumericalParameters)::Tuple{SingleParticleFunction,SingleParticleFunction,Array{Float64,2}}

   nflav=length(pSimulation.gf_flav)
   Nmax=pModel.Nmax
   beta = pSimulation.beta

   gfdiagnorm                    = zeros(Float64,nflav)                   # check normalization of GF
   gf_w::SingleParticleFunction  = zeros(nflav,nflav,length(pFreq.wf))
   gf_iw::SingleParticleFunction = zeros(nflav,nflav,length(pFreq.iwf))

   #evalContributions::Array{Array{Float64,1},1} = []

   E0 = minimum( first.(evallist) )
   evalContributions = Array{Float64}(undef, (length(evallist), 4) )

   # prefill evalContribution array
   for n1=1:length(evallist)
      E1    = evallist[NSperm[n1]][1] -E0 
      N1    = round(Int64,evallist[NSperm[n1]][2])
      s1    = round(Int64,evallist[NSperm[n1]][3])
      S1    = spinConfig(s1,N1,Nmax)
      #push!( evalContributions, [N1,S1,E1,0.0] )
      evalContributions[n1, :] .= [N1,S1,E1,0.0]
   end


   # now we loop over the flavors
   # GF is defined as in the pdf docu:
   # <n1| c_a |n2><n2| cdag_b |n1>
   for a=1:length(pSimulation.gf_flav)
      for b=a:length(pSimulation.gf_flav) # only upper triangular part

         # first count how many transitions we have
         nopairs = 0
         for i=1:length(possibleTransitions[2*a-1])
            n2set = intersect( possibleTransitions[2*a-1][i],possibleTransitions[2*b-1][i] )
            Ei = evallist[NSperm[i]][1] -E0
            for j in n2set
               Ej = evallist[NSperm[j]][1] -E0
               if exp(-beta*Ei)+exp(-beta*Ej) > pNumerics.cutoff
                  nopairs += 1
               end
            end
         end

         nmpairs = zeros(Int64,nopairs,2)
         ii = 1
         for i=1:length(possibleTransitions[2*a-1])
            n2set = intersect( possibleTransitions[2*a-1][i],possibleTransitions[2*b-1][i] )
            Ei = evallist[NSperm[i]][1] -E0
            for j in n2set
               Ej = evallist[NSperm[j]][1] -E0
               if exp(-beta*Ei)+exp(-beta*Ej) > pNumerics.cutoff
                  nmpairs[ii,:] = [i,j]
                  ii += 1
               end
            end
         end

         for i=1:nopairs
            n = nmpairs[i,1]
            m = nmpairs[i,2]

            E1    = evallist[NSperm[n]][1] -E0     # shift by E0 to avoid overflow
            N1    = round(Int64,evallist[NSperm[n]][2])
            s1    = round(Int64,evallist[NSperm[n]][3])
            S1    = spinConfig(s1,N1,Nmax)
   
            E2    = evallist[NSperm[m]][1] -E0
            N2    = round(Int64,evallist[NSperm[m]][2])   
            s2    = round(Int64,evallist[NSperm[m]][3])
            S2    = spinConfig(s2,N2,Nmax)

            # It's faster to calculate exp than to store/read it

            ovrlp = (exp(-beta*E1)+exp(-beta*E2))*overlaps[2*a-0,n,m]*overlaps[2*b-1,m,n]            # <n1|c_a|n2> * <n2|cdag_b|n1>
            gf_w[a,b,:]  += ovrlp ./ ( pFreq.wf     .+ (pNumerics.delta*im + E1 - E2) )
            gf_iw[a,b,:] += ovrlp ./ ( im*pFreq.iwf .+ (                   + E1 - E2) )
            evalContributions[n,4] += abs(ovrlp)
            evalContributions[m,4] += abs(ovrlp)
            if a==b
               gfdiagnorm[a] += real(ovrlp)
            end

         end # n,m loop
      end # b loop
   end # a loop
    println("\rdone!")
   
   # copy the offdiagonals (not true copy in julia, just referencing but this is ok)
   for m1=1:nflav
      for m2=m1+1:nflav
         gf_w[m2,m1,:] = gf_w[m1,m2,:]
         gf_iw[m2,m1,:] = gf_iw[m1,m2,:]
      end
   end
   
   #normalize
   gfdiagnorm ./= getZ(evallist,beta)
   gf_w ./= getZ(evallist,beta)
   gf_iw ./= getZ(evallist,beta)

   # normalize the Matsubara GF tail so that the Dyson equation works for the Selfenergy when the normalization is slightly too small
#   for m1=1:nflav
#      gf_iw[m1,m1,:] ./= gfdiagnorm[m1]
#   end

   writeGFnorm(gfdiagnorm, maximum( first.(evallist) )-E0)

   return gf_w, gf_iw, evalContributions
end 
#########################################################################################
