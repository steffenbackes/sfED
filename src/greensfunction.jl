const SingleParticleFunction = Array{Complex{Float32},3}
const TwoParticleFunction = Array{Complex{Float32},3}  # single orbital for now, depends on 3 frequencies


########################################################################
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
    getGF(evallist, eveclist, allstates, pModel, pSimulation, pFreq, pNumerics)

Evaluation of the Lehmann representation for the interacting finite-temperature Green's function
We sum over all Eigenstates `evallist` with electron number N and all other eigenstates n2 with electron number N-1
Since we have spin degeneracy, we only calculate the up GF, so we annihilate one up spin, so we also only sum over all S-1 states for given S(n1)
"""
function getGF(transitions::Array{Transitions,1},   # two for each flavor
               Z::Float64,
               pSimulation::SimulationParameters,
               pFreq::FrequencyMeshes,
               pNumerics::NumericalParameters)::Tuple{SingleParticleFunction,SingleParticleFunction} #,Array{Float64,2}}

   Nmax=length(transitions[1].transitions)-1
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

   # now we loop over the remaining transitions:  # GF is defined as in the pdf docu:  # <n1| c_a |n2><n2| cdag_b |n1>
   for b=1:length(pSimulation.gf_flav)          # flavor loop
      # loop over all transitions c^dag_b|1>
      for n1=0:Nmax
         for s1=1:noSpinConfig(n1,Nmax)
            for trans2cdagb1 in transitions[2*b-1].transitions[n1+1][s1]
               overlapb = trans2cdagb1.overlap
               ifrom,ito = trans2cdagb1.iFromTo[1:2]
               Efrom,Eto = trans2cdagb1.EvalFromTo[1:2]
               n2 = n1+1
               s2 = trans2cdagb1.sFromTo[2]

               for a=1:b # only upper triangular part
                  # the transition <1|c_a|2> is completely fixed by knowing <2|c^dag_b|1>, just look up the overlap if there exists one
                  for it in findall(x->(x.iFromTo==[ito,ifrom]), transitions[2*a-0].transitions[n2+1][s2] )

                     ovrlp = (exp(-beta*Efrom)+exp(-beta*Eto))*transitions[2*a-0].transitions[n2+1][s2][it].overlap*overlapb      # <n1|c_a|n2> * <n2|cdag_b|n1>

                     if abs(ovrlp)>pNumerics.cutoff
                        gf_w[a,b,:]  += ovrlp ./ ( pFreq.wf     .+ (pNumerics.delta*im + Efrom - Eto) )
                        gf_iw[a,b,:] += ovrlp ./ ( im*pFreq.iwf .+ (                   + Efrom - Eto) )
                       # evalContributions[n,4] += abs(ovrlp)
                       # evalContributions[m,4] += abs(ovrlp)
                       a==b && (gfdiagnorm[a] += real(ovrlp))
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
    getFuckingLarge2partTerm(w1,w2,w3,Em,En,Eo,Ep)
Evaluate the Lehmann term for the 2part Green's function
"""
function getFuckingLarge2partTerm(w1::Float32,w2::Float32,w3::Float32,
                                  Em::Eigenvalue,En::Eigenvalue,Eo::Eigenvalue,Ep::Eigenvalue,
                                  expEop::Float32,expEnp::Float32,expEmp::Float32)::Complex{Float32}

   return (  ( expEop/(im*w3+Eo-Ep) + expEnp/(im*w2+im*w3+En-Ep) )/(im*w2+En-Eo)  
            -( expEop/(im*w3+Eo-Ep) - expEmp/(im*w1+(im*w2+im*w3+Em-Ep)) )/(im*w1+(im*w2+Em-Eo))
          )/( im*w1+(Em-En) )
end

#####################################################

"""
    getGF2partTerm(transitions,pFreq, pNumerics)
Evaluate one of the six terms for the GF, the only difference is where to read out
the overlaps and the permutation of the frequencies
"""
function getGF2partTerm(transitionNtoM::Transitions,
                        transitionOtoN::Transitions,
                        transitionPtoO::Transitions,
                        transitionMtoP::Transitions,
                        pFreq::FrequencyMeshes, nw::Int64,
                        wperm::Array{Int64,1},
                        pNumerics::NumericalParameters)::Tuple{TwoParticleFunction, Float32 }

   gfnorm::Float32 = 0.0
   gfterm = zeros(Complex{Float32},nw,nw,nw)
   Nmax=length(transitionMtoP.transitions)-1

   for nm=0:Nmax
      for sm=1:noSpinConfig(nm,Nmax)

         for transMtoP in transitionMtoP.transitions[nm+1][sm]
            mstate = transMtoP.iFromTo[1]
            Em,Ep = transMtoP.EvalFromTo[1:2]
            expEmp = sum(transMtoP.ExpFromTo[1:2])
            np = transMtoP.nFromTo[2]
            sp = transMtoP.sFromTo[2]
            pstate = transMtoP.iFromTo[2]
      
            # now pick out the <o|c|p> transitions that have p=pstate
            itransPO = findall(x->(x.iFromTo[1]==pstate), transitionPtoO.transitions[np+1][sp] )
            for transPtoO in transitionPtoO.transitions[np+1][sp][itransPO]
               Eo = transPtoO.EvalFromTo[2]
               expEop = sum(transPtoO.ExpFromTo[1:2])
               no = transPtoO.nFromTo[2]
               so = transPtoO.sFromTo[2]
               ostate = transPtoO.iFromTo[2]
      
               # now pick out the <n|c|o> transitions that have o=ostate
               itransON = findall(x->(x.iFromTo[1]==ostate), transitionOtoN.transitions[no+1][so] )
               for transOtoN in transitionOtoN.transitions[no+1][so][itransON]
                  En = transOtoN.EvalFromTo[2]
                  expEnp = transOtoN.ExpFromTo[2] - transPtoO.ExpFromTo[1]
                  nn = transOtoN.nFromTo[2]
                  sn = transOtoN.sFromTo[2]
                  nstate = transOtoN.iFromTo[2]
      
                  # Apply exp cutoff
                  if expEop + expEnp + expEmp > pNumerics.cutoff
      
                     # now pick out the <m|c|n> transitions that have n=nstate AND m=mstate since we need to go back
                     itransNM = findall(x->(x.iFromTo==[nstate,mstate]), transitionNtoM.transitions[nn+1][sn] )
                     for transNtoM in transitionNtoM.transitions[nn+1][sn][itransNM]
      
                        overlap = transMtoP.overlap * transPtoO.overlap * transOtoN.overlap * transNtoM.overlap
                        gfnorm += real(overlap)
                        for n1=1:nw
                           for n2=1:nw
                              for n3=1:nw
                                 w = [ pFreq.iwf[n1], pFreq.iwf[n2], pFreq.iwf[n3] ]
                                 gfterm[n1,n2,n3] += overlap*getFuckingLarge2partTerm(w[wperm[1]],w[wperm[2]],w[wperm[3]],
                                                                                      Em,En,Eo,Ep,
                                                                                      expEop,expEnp,expEmp)
                              end # n3
                           end # n2
                        end # n1
      
                     end # n->m
                  end # exp cutoff
               end # o->n
            end # p->o
         end # m->p

      end # sm
   end # nm
   return gfterm, gfnorm
end

###############################################################################################
"""
    getGF2part(evallist, eveclist, allstates, pModel, pSimulation, pFreq, pNumerics)

Evaluation of the Lehmann representation for the interacting finite-temperature 2-particle Green's function
We sum over all Eigenstates `evallist` with electron number N and all other eigenstates connected by c and c^dagger
We choose the definition G^(2)_up,dn = <T c^dag_up(t1) c_up(t2) c^dag_dn(t3) c_dn(0) >
"""
function getGF2part(transitions::Array{Transitions,1},   # two for each flavor
                    Z::Float64,
                    pFreq::FrequencyMeshes,nw::Int64,
                    pNumerics::NumericalParameters)::TwoParticleFunction        #Tuple{TwoParticleFunction,Array{Float64,2}}

   # we restrict this calculation to one orbital !
   println("Generating the two-particle GF on ",nw,"^3 frequencies...")
   #evalContributions = Array{Float64}(undef, (length(evallist), 4) ) REIMPLEMENT THIS ASAP
   ## prefill the evalContributions array
   #for m=1:length(evallist)
   #   Em    = evallist[NSperm[m]][1] -E0 
   #   Nm    = round(Int64,evallist[NSperm[m]][2])
   #   sm    = round(Int64,evallist[NSperm[m]][3])
   #   Sm    = spinConfig(sm,Nm,Nmax)
   #   evalContributions[m,:] = [Nm,Sm,Em,0.0]
   #end

   gf,gfnorm = .-getGF2partTerm(transitions[1],transitions[2],transitions[3],transitions[4],pFreq,nw,[1,2,3],pNumerics)
               .+getGF2partTerm(transitions[2],transitions[1],transitions[3],transitions[4],pFreq,nw,[2,1,3],pNumerics)
               .-getGF2partTerm(transitions[2],transitions[3],transitions[1],transitions[4],pFreq,nw,[2,3,1],pNumerics)
               .+getGF2partTerm(transitions[3],transitions[2],transitions[1],transitions[4],pFreq,nw,[3,2,1],pNumerics)
               .+getGF2partTerm(transitions[1],transitions[3],transitions[2],transitions[4],pFreq,nw,[1,3,2],pNumerics)
               .-getGF2partTerm(transitions[3],transitions[1],transitions[2],transitions[4],pFreq,nw,[3,1,2],pNumerics)

   println("\rdone!")
   gf ./= Z
   gfnorm /= Z
   @printf("2-particle Green's function normalized to %.3f \n",gfnorm)
   return gf
end 

#######################################################################

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
