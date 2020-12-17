"""
    getG0(eps, tmatrix, w)

Construct the noninteracting Green's function from onsite energy `eps`,
hoppings between sites `i` and `j`, `tmatrix` and Matsubara grid `w`.
"""
function getG0(        eps::Array{Float64,1},
                   tmatrix::Array{Float64,2},  
               pSimulation::SimulationParameters,
                         w::FrequencyMeshCplx )::SingleParticleFunction
   nflav=length(pSimulation.gf_flav)
   nw = length(w)
   gf0::SingleParticleFunction = zeros(nflav,nflav,nw)

   # Since the flavors can be up/dn, we need to blow up the tmatrix, which is only defined for orbitals, to spin/orbital space
   norb = size(tmatrix)[1]
   tmatspin = zeros(2*norb,2*norb)
   for m1=0:norb-1
      for m2=0:norb-1
         for s=0:1
            tmatspin[2*m1+s+1,2*m2+s+1] = tmatrix[m1+1,m2+1]
         end
      end
   end
   for n=1:nw
      gf0[:,:,n] = inv( I*(w[n] + pSimulation.mu) - Diagonal(eps[pSimulation.gf_flav]) - tmatspin[pSimulation.gf_flav,pSimulation.gf_flav] ) 
   end
   return gf0
end
########################################################################


"""
    getCmatrix(anni::Int64,subspace1, subspace2)

Generate the matrix for the annihilation operator acting on the subspace with N,S quantum number
we act `subspace1` * c * `subspace2`
`subspace 2` has all basis states with N,S
`subspace 1` has all basis states with N-1,S-dS
now just determine the transformation when acting the annihilation operator `anni` on it
`anni` is the orbital/spin index
"""
function getCmatrix(     anni::Int64, 
                       subspace1::Array{Fockstate,1}, 
                       subspace2::Array{Fockstate,1})::CAmatrix
   dim1 = length(subspace1)
   dim2 = length(subspace2)
   indexI = Int64[]
   indexJ = Int64[]
   val = Int64[]
   # now annihilate the particle in all basis states in subspace2 and find the corresponding state in subspace1
   for j=1:dim2
      if subspace2[j][anni] == 1   # if there is one particle to annihilate...
         state = copy( subspace2[j] )
         c1sgn = getCsign(anni,state) # count the particles before anni
         state[anni] = 0

         #println("Find ",state," in list of ",subspace1)
         # now find this state in the subspace1
         i = findfirst(map(x-> all(x .== state), subspace1))
         push!(indexI,i)
         push!(indexJ,j)
         push!(val,c1sgn)
      end
   end # j
   return sparse(indexI,indexJ,val, dim1,dim2)
end

"""
    getCdagmatrix(crea,subspace1,subspace2)

Generate the matrix for the creation operator acting on the subspace with N-1,S-1 quantum numbers
we act `subspace1` * cdag * `subspace2`
`subspace 2` has all basis states with N-1,S-1 because we have previously annihilated one up electron
`subspace 1` has all basis states with N,S
now just determine the transformation when acting the creation operator `crea` on it
`crea` is the orbital/spin index
"""
function getCdagmatrix(     crea::Int64, 
                       subspace1::Array{Fockstate,1}, 
                       subspace2::Array{Fockstate,1})::CAmatrix
   dim1 = length(subspace1)
   dim2 = length(subspace2)
   indexI = Int64[]
   indexJ = Int64[]
   val = Int64[]

   # now create the particle in all basis states and find the corresponding state in the N,S space
   for j=1:dim2
      if subspace2[j][crea] == 0   # if there is space to create one particle...
         state = copy( subspace2[j] )
         c1sgn = getCsign(crea,state) # count the particles before crea
         state[crea] = 1

         #println("Find ",state," in list of ",subspace1)
         # now find this state in the N,S subspace
         i = findfirst(map(x-> all(x .== state), subspace1))
         push!(indexI,i)
         push!(indexJ,j)
         push!(val,c1sgn)
      end
   end # j
   return sparse(indexI,indexJ,val, dim1,dim2)
end
#######################################################################


"""
    getGF(evallist, eveclist, allstates, pModel, pSimulation, pFreq, pNumerics)

Evaluation of the Lehmann representation for the interacting finite-temperature Green's function
We sum over all Eigenstates `evallist` with electron number N and all other eigenstates n2 with electron number N-1
Since we have spin degeneracy, we only calculate the up GF, so we annihilate one up spin, so we also only sum over all S-1 states for given S(n1)
"""
function getGF(transitions::Array{Transitions,1},
               Z::Float64,
               pSimulation::SimulationParameters,
               pFreq::FrequencyMeshes,
               pNumerics::NumericalParameters)::Tuple{SingleParticleFunction,SingleParticleFunction} #,Array{Float64,2}}

   nflav=length(pSimulation.gf_flav)
   beta = pSimulation.beta

   gfdiagnorm                    = zeros(Float64,nflav)                   # check normalization of GF
   gf_w::SingleParticleFunction  = zeros(nflav,nflav,length(pFreq.wf))
   gf_iw::SingleParticleFunction = zeros(nflav,nflav,length(pFreq.iwf))
   #evalContributions = Array{Float64}(undef, (length(evallist), 4) )

   # prefill evalContribution array     # PUT ON HOLD FOR NOW, BUT SHOULD BE REIMPLEMENTED ASAP!
   #for n1=1:length(evallist)
   #   E1    = evallist[NSperm[n1]][1] -E0 
   #   N1    = round(Int64,evallist[NSperm[n1]][2])
   #   s1    = round(Int64,evallist[NSperm[n1]][3])
   #   S1    = spinConfig(s1,N1,Nmax)
   #   #push!( evalContributions, [N1,S1,E1,0.0] )
   #   evalContributions[n1, :] .= [N1,S1,E1,0.0]
   #end

   # now we loop over the remaining transitions
   # GF is defined as in the pdf docu:
   # <n1| c_a |n2><n2| cdag_b |n1>
   for b=1:length(pSimulation.gf_flav)
      # loop over all transitions c^dag_b|1>
      transCdagB = transitions[2*b-1]

      for n1=0:transCdagB.Nmax
         for s1=1:transCdagB.nstatesS[n1+1]
            for i=1:transCdagB.dimNS[n1+1,s1]
               t2cdagb1 = transCdagB[n1,s1,i]
               ifrom = Int64(t2cdagb1[1])
               ito   = Int64(t2cdagb1[2])
               Efrom = t2cdagb1[3]
               Eto   = t2cdagb1[4]
               isfrom = Int64(t2cdagb1[5])
               isto   = Int64(t2cdagb1[6])
               ovrlpB = t2cdagb1[7]
               nto = n1+1

               for a=1:b # only upper triangular part
                  transCA = transitions[2*a-0]  # <1|c_a|2> # the transition is already fixed by <n2| |n1>

                  ovrlp = (exp(-beta*Efrom)+exp(-beta*Eto))*transCA[nto,isto,ito,ifrom]*ovrlpB      # <n1|c_a|n2> * <n2|cdag_b|n1>

                  if abs(ovrlp)>pNumerics.cutoff
                     gf_w[a,b,:]  += ovrlp ./ ( pFreq.wf     .+ (pNumerics.delta*im + Efrom - Eto) )
                     gf_iw[a,b,:] += ovrlp ./ ( im*pFreq.iwf .+ (                   + Efrom - Eto) )
                    # evalContributions[n,4] += abs(ovrlp)
                    # evalContributions[m,4] += abs(ovrlp)
                     if a==b
                        gfdiagnorm[a] += real(ovrlp)
                     end
                  end
               end # a loop
            end # i transition from <2|c^dag_b|1>
         end # s1
      end # n1 loop
   end # b loop
    println("\rdone!")
   
   # copy the offdiagonals (not true copy in julia, just referencing but this is ok)
   for m1=1:nflav
      for m2=m1+1:nflav
         gf_w[m2,m1,:] = gf_w[m1,m2,:]
         gf_iw[m2,m1,:] = gf_iw[m1,m2,:]
      end
   end
   
   #normalize
   gfdiagnorm ./= Z
   gf_w ./= Z
   gf_iw ./= Z

   # normalize the Matsubara GF tail so that the Dyson equation works for the Selfenergy when the normalization is slightly too small
   # Not needed since we diagonalize exactly
#   for m1=1:nflav
#      gf_iw[m1,m1,:] ./= gfdiagnorm[m1]
#   end

   writeGFnorm(gfdiagnorm)

   return gf_w, gf_iw #, evalContributions
end 
#########################################################################################
"""
    getFuckingLarge2partTerm(w1,w2,w3,Em,En,Eo,Ep,beta)
Evaluate the Lehmann term for the 2part Green's function
"""
function getFuckingLarge2partTerm(w1::FrequencyMesh,w2::Float32,w3::Float32,Em,En,Eo,Ep,beta)::Array{Complex{Float64},1}
   expEop = exp(-beta*Eo) + exp(-beta*Ep)
   expEnp = exp(-beta*En) - exp(-beta*Ep)
   expEmp = exp(-beta*Em) + exp(-beta*Ep)

   return (  ( expEop/(im*w3+Eo-Ep) + expEnp/(im*w2+im*w3+En-Ep) )/(im*w2+En-Eo)  
           .-( expEop/(im*w3+Eo-Ep) .- expEmp./(im.*w1.+(im*w2+im*w3+Em-Ep)) )./(im.*w1.+(im*w2+Em-Eo))
          )./( im.*w1.+(Em-En) )
end

###############################################################################################
"""
    getGF2part(evallist, eveclist, allstates, pModel, pSimulation, pFreq, pNumerics)

Evaluation of the Lehmann representation for the interacting finite-temperature 2-particle Green's function
We sum over all Eigenstates `evallist` with electron number N and all other eigenstates connected by c and c^dagger
We choose the definition G^(2)_up,dn = <T c^dag_up(t1) c_up(t2) c^dag_dn(t3) c_dn(0) >
"""
function getGF2part( evallist::Array{Array{Eigenvalue,1},1},
                eveclist::Array{Eigenvector,1},
               allstates::NSstates,
               NSperm::Array{Int64,1},
               pSimulation::SimulationParameters,
               pFreq::FrequencyMeshes,
               pNumerics::NumericalParameters)::Tuple{TwoParticleFunction,Array{Float64,2}}

   # we restrict this calculation to one orbital !
   #norb=1
   println("ERROR: getGF2part NOT IMPLEMENTED")
   exit()
   Nmax=-1
   beta = pSimulation.beta
   nw = 5
   println("Generating the two-particle GF on ",nw,"^3 frequencies...")

   # arrays to check when we need to recalculate the cdagger cmatrices
   NSpm   = Int64[-1,-1,-1,-1]
   NSop13   = Int64[-1,-1,-1,-1]
   NSop46   = Int64[-1,-1,-1,-1]
   NSop25   = Int64[-1,-1,-1,-1]

   NSno1   = Int64[-1,-1,-1,-1]
   NSno3   = Int64[-1,-1,-1,-1]
   NSno4   = Int64[-1,-1,-1,-1]
   NSno6   = Int64[-1,-1,-1,-1]
   NSno2   = Int64[-1,-1,-1,-1]
   NSno5   = Int64[-1,-1,-1,-1]

   NSmn1   = Int64[-1,-1,-1,-1]
   NSmn3   = Int64[-1,-1,-1,-1]
   NSmn4   = Int64[-1,-1,-1,-1]
   NSmn6   = Int64[-1,-1,-1,-1]
   NSmn2   = Int64[-1,-1,-1,-1]
   NSmn5   = Int64[-1,-1,-1,-1]

   # Term 1 creation/annihilation operators
   cmat_dn_pm    = CAmatrix[]              # cmat_dn_pm is used by all, only define here

   cdmat13_dn_op   = CAmatrix[]
   cmat1_up_no    = CAmatrix[]
   cdmat1_up_mn   = CAmatrix[]
   # Term 3 creation/annihilation operators shares p->o with Term 1
   cdmat3_up_no   = CAmatrix[]
   cmat3_up_mn    = CAmatrix[]

   # Term 4 creation/annihilation operators
   cdmat46_up_op   = CAmatrix[]
   cdmat4_dn_no    = CAmatrix[]
   cmat4_up_mn    = CAmatrix[]
   # Term 6 creation/annihilation operators shares p_o with Term 4
   cmat6_up_no    = CAmatrix[]
   cdmat6_dn_mn    = CAmatrix[]

   # Term 2 creation/annihilation operators
   cmat25_up_op   = CAmatrix[]
   cdmat2_dn_no    = CAmatrix[]
   cdmat2_up_mn    = CAmatrix[]
   # Term 5 creation/annihilation operators shares p_o with Term 2
   cdmat5_up_no    = CAmatrix[]
   cdmat5_dn_mn    = CAmatrix[]

   gf::TwoParticleFunction  = zeros(nw,nw,nw)
   gfnorm = 0.0

   actuallyContributedSth = 0

   E0 = minimum( first.(evallist) )

   evalContributions = Array{Float64}(undef, (length(evallist), 4) )
   # prefill the evalContributions array
   for m=1:length(evallist)
      Em    = evallist[NSperm[m]][1] -E0 
      Nm    = round(Int64,evallist[NSperm[m]][2])
      sm    = round(Int64,evallist[NSperm[m]][3])
      Sm    = spinConfig(sm,Nm,Nmax)
      evalContributions[m,:] = [Nm,Sm,Em,0.0]
   end

   # Now calculate the 2particle GF
   for m=1:length(evallist)
      if m%(max(1,round(Int64,length(evallist)/100.0))) == 0
            print("\r"*lpad(round(Int64,m*100.0/length(evallist)),4)*"%")
      end

      Em    = evallist[NSperm[m]][1] -E0     # shift by E0 to avoid overflow
      Nm    = round(Int64,evallist[NSperm[m]][2])
      sm    = round(Int64,evallist[NSperm[m]][3])
      evecm = eveclist[NSperm[m]]
      Sm    = spinConfig(sm,Nm,Nmax)
      expFm = exp(-beta*Em) # exponential Boltzmann factor

      for n=1:length(evallist)
         En    = evallist[NSperm[n]][1] -E0
         Nn    = round(Int64,evallist[NSperm[n]][2])   
         sn    = round(Int64,evallist[NSperm[n]][3])
         evecn = eveclist[NSperm[n]]
         Sn    = spinConfig(sn,Nn,Nmax)
         expFn = exp(-beta*En)

         for o=1:length(evallist)
            Eo    = evallist[NSperm[o]][1] -E0
            No    = round(Int64,evallist[NSperm[o]][2])   
            so    = round(Int64,evallist[NSperm[o]][3])
            eveco = eveclist[NSperm[o]]
            So    = spinConfig(so,No,Nmax)
            expFo = exp(-beta*Eo)

            for p=1:length(evallist)
               Ep    = evallist[NSperm[p]][1] -E0
               Np    = round(Int64,evallist[NSperm[p]][2])   
               sp    = round(Int64,evallist[NSperm[p]][3])
               evecp = eveclist[NSperm[p]]
               Sp    = spinConfig(sp,Np,Nmax)

               expFop = expFo+exp(-beta*Ep)
               expFmp = expFm+exp(-beta*Ep)
               expFnp = expFn-exp(-beta*Ep)
         
               # Exclude transitions too high in energy where all terms are zero and check for right m->p transition
               if ( expFop+expFmp+expFnp > pNumerics.cutoff && Np==Nm-1 && Sp==Sm+1)
         
                  # Check for new N,S combinations at which creation/annihilation matrices we need to update, m->p is used by all terms
                  if ( NSpm!=[Np,Sp, Nm,Sm] )
                     NSpm = [Np,Sp, Nm,Sm]
                     cmat_dn_pm = getCmatrix(2, allstates[Np+1][sp], allstates[Nm+1][sm]) 
                  end
                     
                  ########################################################################################################
                  # Term 1 + 3 as defined in the pdf
                  if ( No==Np+1 && So==Sp-1) # check p->o transition which is the same for 1+3
                     # recalculate matrix used by Term 1 + 3
                     if ( NSop13!=[No,So, Np,Sp] )
                        NSop13=[No,So, Np,Sp]
                        cdmat13_dn_op = getCdagmatrix(2, allstates[No+1][so], allstates[Np+1][sp]) 
                     end
                     # This overlap is the same for 1+3
                     overlapopm = dot( eveco, cdmat13_dn_op*evecp )*dot( evecp, cmat_dn_pm*evecm )

                     if abs(overlapopm)>pNumerics.cutoff

                        # Term 1 ###############################################################################
                        if ( Nn==No-1 && Sn==So-1) # check o->n, n->m transition is automatically fulfilled then
                           # recalculate matrix
                           if ( NSno1!=[Nn,Sn, No,So] )
                              NSno1=[Nn,Sn, No,So]
                              cmat1_up_no = getCmatrix(1, allstates[Nn+1][sn], allstates[No+1][so]) 
                           end
                           if ( NSmn1!=[Nm,Sm, Nn,Sn] )
                              NSmn1=[Nm,Sm, Nn,Sn]
                              cdmat1_up_mn = getCdagmatrix(1, allstates[Nm+1][sm], allstates[Nn+1][sn]) 
                           end
                           # now we have everything
                           overlap = dot( evecm, cdmat1_up_mn*evecn )*dot( evecn, cmat1_up_no*eveco )*overlapopm
                           if abs(overlap)>pNumerics.cutoff
                              actuallyContributedSth += 1
                              evalContributions[m,4] += abs(overlap)
                              evalContributions[n,4] += abs(overlap)
                              evalContributions[o,4] += abs(overlap)
                              evalContributions[p,4] += abs(overlap)
                              gfnorm += -real(overlap)
                              for w2=1:nw    # Here loop over all frequencies, the first argument is vectorized
                                 for w3=1:nw
                                    gf[:,w2,w3] += -getFuckingLarge2partTerm(pFreq.iwf[1:nw],pFreq.iwf[w2],pFreq.iwf[w3],Em,En,Eo,Ep,beta) .*overlap
                                 end
                              end
                           end # overlap cutoff Term 1
                        end #p->n Term 1

                        # Term 3 ###############################################################################
                        if ( Nn==No+1 && Sn==So+1) # check o->n, n->m transition is automatically fulfilled then
                           # recalculate matrix
                           if ( NSno3!=[Nn,Sn, No,So] )
                              NSno3=[Nn,Sn, No,So]
                              cdmat3_up_no = getCdagmatrix(1, allstates[Nn+1][sn], allstates[No+1][so])
                           end
                           if ( NSmn3!=[Nm,Sm, Nn,Sn] )
                              NSmn3!=[Nm,Sm, Nn,Sn]
                              cmat3_up_mn = getCmatrix(1, allstates[Nm+1][sm], allstates[Nn+1][sn])
                           end
                           # now we have everything
                           overlap = dot( evecm, cmat3_up_mn*evecn )*dot( evecn, cdmat3_up_no*eveco )*overlapopm
                           if abs(overlap)>pNumerics.cutoff
                              actuallyContributedSth += 1
                              evalContributions[m,4] += abs(overlap)
                              evalContributions[n,4] += abs(overlap)
                              evalContributions[o,4] += abs(overlap)
                              evalContributions[p,4] += abs(overlap)
                              gfnorm += real(overlap)
                              for w1=1:nw    # Here loop over all frequencies, the first argument is vectorized
                                 for w3=1:nw
                                    gf[w1,:,w3] += getFuckingLarge2partTerm(pFreq.iwf[1:nw],pFreq.iwf[w1],pFreq.iwf[w3],Em,En,Eo,Ep,beta) .*overlap
                                 end
                              end
                           end # overlap cutoff Term 3
                        end #p->n Term 3

                     end # first overlap cutoff for Term 1+3 
                  end # p->o Transition for 1+3
                  ########################################################################################################

                  ########################################################################################################
                  # Term 4 + 6 as defined in the pdf
                  if ( No==Np+1 && So==Sp+1) # check p->o transition which is the same for 4+6
                     # recalculate matrix used by Term 4 + 6
                     if ( NSop46!=[No,So, Np,Sp] )
                        NSop46=[No,So, Np,Sp]
                        cdmat46_up_op = getCdagmatrix(1, allstates[No+1][so], allstates[Np+1][sp])
                     end
                     # This overlap is the same for 4+6
                     overlapopm = dot( eveco, cdmat46_up_op*evecp )*dot( evecp, cmat_dn_pm*evecm )

                     if abs(overlapopm)>pNumerics.cutoff

                        # Term 4 ###############################################################################
                        if ( Nn==No+1 && Sn==So-1) # check o->n, n->m transition is automatically fulfilled then
                           # recalculate matrix
                           if ( NSno4!=[Nn,Sn, No,So] )
                              NSno4=[Nn,Sn, No,So]
                              cdmat4_dn_no = getCdagmatrix(2, allstates[Nn+1][sn], allstates[No+1][so])
                           end
                           if ( NSmn4!=[Nm,Sm, Nn,Sn] )
                              NSmn4=[Nm,Sm, Nn,Sn]
                              cmat4_up_mn = getCmatrix(1, allstates[Nm+1][sm], allstates[Nn+1][sn])
                           end
                           # now we have everything
                           overlap = dot( evecm, cmat4_up_mn*evecn )*dot( evecn, cdmat4_dn_no*eveco )*overlapopm
                           if abs(overlap)>pNumerics.cutoff
                              actuallyContributedSth += 1
                              evalContributions[m,4] += abs(overlap)
                              evalContributions[n,4] += abs(overlap)
                              evalContributions[o,4] += abs(overlap)
                              evalContributions[p,4] += abs(overlap)
                              gfnorm += -real(overlap)
                              for w3=1:nw    # Here loop over all frequencies, the first argument is vectorized
                                 for w1=1:nw
                                    gf[w1,:,w3] += -getFuckingLarge2partTerm(pFreq.iwf[1:nw],pFreq.iwf[w3],pFreq.iwf[w1],Em,En,Eo,Ep,beta) .*overlap
                                 end
                              end
                           end # overlap cutoff Term 4
                        end #p->n Term 1

                        # Term 6 ###############################################################################
                        if ( Nn==No-1 && Sn==So-1) # check o->n, n->m transition is automatically fulfilled then
                           # recalculate matrix
                           if ( NSno6!=[Nn,Sn, No,So] )
                              NSno6=[Nn,Sn, No,So]
                              cmat6_up_no = getCmatrix(1, allstates[Nn+1][sn], allstates[No+1][so])
                           end
                           if ( NSmn6!=[Nm,Sm, Nn,Sn] )
                              NSmn6!=[Nm,Sm, Nn,Sn]
                              cdmat6_dn_mn = getCdagmatrix(2, allstates[Nm+1][sm], allstates[Nn+1][sn])
                           end
                           # now we have everything
                           overlap = dot( evecm, cdmat6_dn_mn*evecn )*dot( evecn, cmat6_up_no*eveco )*overlapopm
                           if abs(overlap)>pNumerics.cutoff
                              actuallyContributedSth += 1
                              evalContributions[m,4] += abs(overlap)
                              evalContributions[n,4] += abs(overlap)
                              evalContributions[o,4] += abs(overlap)
                              evalContributions[p,4] += abs(overlap)
                              gfnorm += real(overlap)
                              for w2=1:nw    # Here loop over all frequencies, the first argument is vectorized
                                 for w1=1:nw
                                    gf[w1,w2,:] += getFuckingLarge2partTerm(pFreq.iwf[1:nw],pFreq.iwf[w2],pFreq.iwf[w1],Em,En,Eo,Ep,beta) .*overlap
                                 end
                              end
                           end # overlap cutoff Term 6
                        end #p->n Term 6

                     end # first overlap cutoff for Term 4+6
                  end # p->o Transition for 4+6
                  ########################################################################################################

                  ########################################################################################################
                  # Term 2 + 5 as defined in the pdf
                  if ( No==Np-1 && So==Sp-1) # check p->o transition which is the same for 2+5
                     # recalculate matrix used by Term 2+5
                     if ( NSop25!=[No,So, Np,Sp] )
                        NSop25=[No,So, Np,Sp]
                        cmat25_up_op = getCmatrix(1, allstates[No+1][so], allstates[Np+1][sp])
                     end
                     # This overlap is the same for 2+5
                     overlapopm = dot( eveco, cmat25_up_op*evecp )*dot( evecp, cmat_dn_pm*evecm )

                     if abs(overlapopm)>pNumerics.cutoff

                        # Term 2 ###############################################################################
                        if ( Nn==No+1 && Sn==So-1) # check o->n, n->m transition is automatically fulfilled then
                           # recalculate matrix
                           if ( NSno2!=[Nn,Sn, No,So] )
                              NSno2=[Nn,Sn, No,So]
                              cdmat2_dn_no = getCdagmatrix(2, allstates[Nn+1][sn], allstates[No+1][so])
                           end
                           if ( NSmn2!=[Nm,Sm, Nn,Sn] )
                              NSmn2=[Nm,Sm, Nn,Sn]
                              cdmat2_up_mn = getCdagmatrix(1, allstates[Nm+1][sm], allstates[Nn+1][sn])
                           end
                           # now we have everything
                           overlap = dot( evecm, cdmat2_up_mn*evecn )*dot( evecn, cdmat2_dn_no*eveco )*overlapopm
                           if abs(overlap)>pNumerics.cutoff
                              actuallyContributedSth += 1
                              evalContributions[m,4] += abs(overlap)
                              evalContributions[n,4] += abs(overlap)
                              evalContributions[o,4] += abs(overlap)
                              evalContributions[p,4] += abs(overlap)
                              gfnorm += real(overlap)
                              for w3=1:nw    # Here loop over all frequencies, the first argument is vectorized
                                 for w2=1:nw
                                    gf[:,w2,w3] += getFuckingLarge2partTerm(pFreq.iwf[1:nw],pFreq.iwf[w3],pFreq.iwf[w2],Em,En,Eo,Ep,beta) .*overlap
                                 end
                              end
                           end # overlap cutoff Term 2
                        end #p->n Term 2

                        # Term 5 ###############################################################################
                        if ( Nn==No+1 && Sn==So+1) # check o->n, n->m transition is automatically fulfilled then
                           # recalculate matrix
                           if ( NSno5!=[Nn,Sn, No,So] )
                              NSno5=[Nn,Sn, No,So]
                              cdmat5_up_no = getCdagmatrix(1, allstates[Nn+1][sn], allstates[No+1][so])
                           end
                           if ( NSmn5!=[Nm,Sm, Nn,Sn] )
                              NSmn5!=[Nm,Sm, Nn,Sn]
                              cdmat5_dn_mn = getCdagmatrix(2, allstates[Nm+1][sm], allstates[Nn+1][sn])
                           end
                           # now we have everything
                           overlap = dot( evecm, cdmat5_dn_mn*evecn )*dot( evecn, cdmat5_up_no*eveco )*overlapopm
                           if abs(overlap)>pNumerics.cutoff
                              actuallyContributedSth += 1
                              evalContributions[m,4] += abs(overlap)
                              evalContributions[n,4] += abs(overlap)
                              evalContributions[o,4] += abs(overlap)
                              evalContributions[p,4] += abs(overlap)
                              gfnorm += -real(overlap)
                              for w1=1:nw    # Here loop over all frequencies, the first argument is vectorized
                                 for w2=1:nw
                                    gf[w1,w2,:] += -getFuckingLarge2partTerm(pFreq.iwf[1:nw],pFreq.iwf[w1],pFreq.iwf[w2],Em,En,Eo,Ep,beta) .*overlap
                                 end
                              end
                           end # overlap cutoff Term 5
                        end #p->n Term 5

                     end # first overlap cutoff for Term 2+5
                  end # p->o Transition for 2+5
                  ########################################################################################################

               end # if exp(-beta*E) > cutoff for all terms and m->p transition for all terms

            end # p loop
         end # o loop
      end # n loop
   end # m loop
    println("\rdone!")
   
   #normalize
   gf ./= getZ(evallist,beta)
   gfnorm /= getZ(evallist,beta)

   @printf("2-particle Green's function normalized to %.3f \n",gfnorm)
   @printf("In the sum over %i^4=%.1e Eigenstates only %i terms contributed ( %.1e %%) \n",
           length(evallist),(length(evallist)*1.0)^4,actuallyContributedSth,actuallyContributedSth*100*(1.0/length(evallist))^4 )

   return gf,evalContributions
end 

################################################################################################################
#""" BUGGY AND SLOWER !
#    getN2Transitions(flavor,eigenspace,fockstates,pNumerics,expCutoff)
#
#Determine 
#"""
#function getN2Transitions(n1::Int64,s1::Int64,flavor::Int64,
#                          eigenspace::Eigenspace,
#                          fockstates::Fockstates,
#                          beta::Float64,
#                          pNumerics::NumericalParameters,
#                          expCutoff::Int64)
#   Nmax = fockstates.Nmax
#   E0 = minimum(eigenspace.evals)
#
#   transitions = Transition[]
#
#   # if flavor is positive, calculate c^dagger, if negative calculate c
#   S1 = eigenspace.SpinNS[n1+1,s1] # spin S in this subspace
#   dS = 2*(abs(flavor)%2)-1
#   n2 = n1+sign(flavor)
#   S2 = S1+dS*sign(flavor)
#
#   n2states = nstatesNS(fockstates,n2,S2) # are there any states in the new subspace?
#
#   cmat = sparse([1],[1],[1])
#   if n2states>0
#      s2 = indexSpinConfig(S2,n2,Nmax)       # this works because we checked this subspace exists
#      if flavor<0
#         cmat = getCmatrix(abs(flavor), fockstates[n2,s2,:], fockstates[n1,s1,:])
#      else
#         cmat = getCdagmatrix(abs(flavor), fockstates[n2,s2,:], fockstates[n1,s1,:])
#      end
#
#      for i=1:eigenspace.dimNS[n1+1,s1] # now loop over the |1> subspace
#         E1,evec1 = eigenspace[n1,s1,i]
#         E1-=E0
#
#         for j=1:eigenspace.dimNS[n2+1,s2] # and loop over the <2| subspace
#            E2,evec2 = eigenspace[n2,s2,j]
#            E2-=E0
#
#            if expCutoff==0 || (exp(-beta*E1)+exp(-beta*E2))>pNumerics.cutoff
#               ovrlp = dot( evec2, cmat * evec1 )
#               if abs(ovrlp) > pNumerics.cutoff
#                  index1 = eigenspace.startNSeval[n1+1,s1,i]
#                  index2 = eigenspace.startNSeval[n2+1,s2,j]
#                  push!( transitions, Transition(index1,index2,E1,E2,ovrlp) )
#               end
#            end # exp cutoff if needed
#         end # j states <2|
#      end # i states |1>
#   end # n2states >0
#
#   return transitions
#end

#########################################################################################

"""
    getPossibleTransitions(eigenspace,fockstates,flavors,pNumerics,expCutoff)

Calculate the overlap elements between all Eigenstates acting on c/c^dagger specified by flavors.
The function returns the overlap array and possibleTransitions for cdagger,c for 
the flavors specified in the flavors array (1=orb1,up, 2=orb1,dn, 3=orb2,up, ....)
"""
function Transitions(flavor::Int64,
                     eigenspace::Eigenspace,
                     fockstates::Fockstates,
                     beta::Float64,
                     pNumerics::NumericalParameters,
                     expCutoff::Int64)
	# if flavor>0 we calculate c^dagger
	# if flavor<0 we calculate c

   Nmax = fockstates.Nmax
   nstates = fockstates.Nstates
   E0 = minimum(eigenspace.evals)

   FromToList = [ ]
   EvalFromToList = [ ]
   isFromToList = [ ]
   overlapList = Complex{Float64}[ ]
   nSmax = size(fockstates.startNS)[2]
   dimNS = zeros(Int64,Nmax+1,nSmax)
   dimMax = 0

   nTransitions = 0
   for n1=0:Nmax
      #if n1%(max(1,round(Int64,Nmax/100.0))) == 0
      #      print("\r"*lpad(round(Int64,n1*100.0/Nmax),4)*"%")
      #end
      
      for s1=1:fockstates.nstatesS[n1+1]
         S1 = eigenspace.SpinNS[n1+1,s1] # spin S in this subspace

         # for given N1,S1 we have the same c,c^dag matrices
         dS = 2*(abs(flavor)%2)-1 # spin increases by 1 for crea if up flavor, etc
 
         # #  <2|cmat|1> ##############################################
         n2 = n1 + sign(flavor)
         S2 = S1 + dS*sign(flavor)
         n2states = nstatesNS(fockstates,n2,S2) # are there any states in the new subspace?

         if n2states>0
            s2 = indexSpinConfig(S2,n2,Nmax)       # this works because we checked this subspace exists
				cmat = sparse([1],[1],[1])
				if flavor>0
	            cmat = getCdagmatrix(abs(flavor), fockstates[n2,s2,:], fockstates[n1,s1,:])
            else
	            cmat = getCmatrix(abs(flavor), fockstates[n2,s2,:], fockstates[n1,s1,:])
            end

            for i=1:eigenspace.dimNS[n1+1,s1] # now loop over the |1> subspace
               Efrom,evecfrom = eigenspace[n1,s1,i]
               Efrom-=E0

               for j=1:eigenspace.dimNS[n2+1,s2] # and loop over the <2| subspace
                  Eto,evecto = eigenspace[n2,s2,j]
                  Eto-=E0

                  if expCutoff==0 || (exp(-beta*Efrom)+exp(-beta*Eto))>pNumerics.cutoff
                     ovrlp = dot( evecto, cmat * evecfrom )
                     if abs(ovrlp) > pNumerics.cutoff
                        ifrom = eigenspace.startNSeval[n1+1,s1,i]
                        ito = eigenspace.startNSeval[n2+1,s2,j]
                        push!( FromToList, [ifrom,ito] )
                        push!( EvalFromToList, [Efrom,Eto] )
                        push!( isFromToList, [s1,s2] )
                        push!( overlapList, ovrlp )

                        dimNS[n1+1,s1] += 1
                        nTransitions += 1
                        if dimNS[n1+1,s1]>dimMax
                           dimMax = dimNS[n1+1,s1]
                        end
                        #println(index1,"->",index2,": N1=",n1,", S1=",S1,", N2=",n2,", S2=",S2)
                     end # overlap cutoff
                  end # exp cutoff if needed
               end # j states <2|
            end # i states |1>
         end # n2states >0

      end # s1
   end # n1
   #println("\rdone!")

   startNS    = zeros(Int64,Nmax+1,nSmax,dimMax)
   FromTo     = Array{Int64}(undef,nTransitions,2)
   isFromTo   = Array{Int64}(undef,nTransitions,2)
   EvalFromTo = Array{Eigenvalue}(undef,nTransitions,2)
   ii=1
   for n1=0:Nmax
      for s1=1:fockstates.nstatesS[n1+1]
         for i=1:dimNS[n1+1,s1]
            startNS[n1+1,s1,i] = ii
            FromTo[ii,:] .= FromToList[ii]
            EvalFromTo[ii,:] .= EvalFromToList[ii]
            isFromTo[ii,:] .= isFromToList[ii]
            ii+=1
         end
      end
   end
   
   return Transitions(Nmax,fockstates.nstatesS,dimNS,startNS,FromTo,EvalFromTo,isFromTo,overlapList)
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
"""
    getSigma(G0, G)

Calculate the Selfenergy from `G0` and `G`
"""
function getSigma(G0::SingleParticleFunction, G::SingleParticleFunction)
   sigma::SingleParticleFunction = zeros(size(G))
   for i=1:size(G)[3]
      sigma[:,:,i]  = inv( G0[:,:,i]) - inv( G[:,:,i])
   end
   return sigma
end
