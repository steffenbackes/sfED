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
	tmatspin[1:norb,1:norb] = tmatrix
	tmatspin[norb+1:end,norb+1:end] = tmatrix
	for n=1:nw
		gf0[:,:,n] = inv( I*(w[n] + pSimulation.mu) - Diagonal(eps[pSimulation.gf_flav]) - tmatspin[pSimulation.gf_flav,pSimulation.gf_flav] ) 
	end
	return gf0
end
########################################################################


"""
    getCmatrix(anni::Int64,subspace2, subspace1)

Generate the matrix for the annihilation operator acting on the subspace with N,S quantum number
we act `subspace2` * c * `subspace1`
`subspace 1` has all basis states with N,S
`subspace 2` has all basis states with N-1,S-1 because we always annihilate one up electron
now just determine the transformation when acting the annihilation operator `anni` on it
`anni` is the orbital/spin index
"""
function getCmatrix(     anni::Int64, 
                    subspace2::Array{Fockstate,1}, 
						  subspace1::Array{Fockstate,1})::CAmatrix
	dim1 = length(subspace1)
	dim2 = length(subspace2)
	indexI = Int64[]
	indexJ = Int64[]
	val = Int64[]

	# now annihilate the particle in all basis states and find the corresponding state in the N-1,S-1 space
	for j=1:dim1
		if subspace1[j][anni] == 1   # if there is a particle to annihilate..
			state = copy( subspace1[j] )
			state[anni] = 0
			a1sgn = (-1)^sum(state[1:anni]) # count the particles before anni

			# now find this state in the smaller subspace
			i = findfirst(map(x-> all(x .== state), subspace2))
			push!(indexI,i)
			push!(indexJ,j)
			push!(val,a1sgn)
		end
	end # j
	return sparse(indexI,indexJ,val, dim2,dim1)
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
			c1sgn = (-1)^sum(state[1:crea]) # count the particles before crea
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
function getGFold( evallist::Array{Array{Eigenvalue,1},1},
                eveclist::Array{Eigenvector,1},
					allstates::NSstates,
					pModel::ModelParameters,
					pSimulation::SimulationParameters,
					pFreq::FrequencyMeshes,
					pNumerics::NumericalParameters)::Tuple{SingleParticleFunction,SingleParticleFunction,Array{Array{Float64,1},1}}

	norb=length(pSimulation.gf_flav)
	Nmax=UInt32(pModel.Nmax)
	beta = pSimulation.beta

	NSvalues   = Int64[-1,-1,-1,-1]
	cmatrix    = CAmatrix[]
	cdagmatrix = CAmatrix[]
	gfdiagnorm                    = zeros(Float64,norb)                   # check normalization of GF
	gf_w::SingleParticleFunction  = zeros(norb,norb,length(pFreq.wf))
	gf_iw::SingleParticleFunction = zeros(norb,norb,length(pFreq.iwf))

	evalContributions::Array{Array{Float64,1},1} = []

	# First create a sortting for all Eigenstates in order of increasing N,S, which makes the GF generation more efficient
	NSperm = getNSperm(evallist)

	E0 = minimum( first.(evallist) )

	for n1=1:length(evallist)
	    if n1%(round(Int64,length(evallist)/10)) == 0
            print("\r"*lpad(round(Int64,n1*100.0/length(evallist)),4)*"%")
		end

		E1    = evallist[NSperm[n1]][1] -E0     # shift by E0 to avoid overflow
		N1    = round(Int64,evallist[NSperm[n1]][2])
		s1    = round(UInt64,evallist[NSperm[n1]][3])
		evec1 = eveclist[NSperm[n1]]
		S1    = spinConfig(s1,N1,Nmax)
		push!( evalContributions, [N1,S1,E1,0.0] )
	
		for n2=1:length(evallist)
			E2    = evallist[NSperm[n2]][1] -E0
			N2    = round(Int64,evallist[NSperm[n2]][2])   
			s2    = round(UInt64,evallist[NSperm[n2]][3])
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
							evalContributions[n1][4] += abs(overlap)
							evalContributions[n2][4] += abs(overlap)
							gfdiagnorm[m1] += real(overlap)
						end
		
						# then upper offdiagonal
						for m2=m1+1:norb
							overlap =  evec1cdagevec2 * dot( evec2, cmatrix[m2]*evec1 ) * expFac # include the Boltzmann terms
							if abs(overlap)>pNumerics.cutoff
								gf_w[m1,m2,:]  += overlap ./ ( pFreq.wf     .+ (pNumerics.delta*im - E1 + E2) )
								gf_iw[m1,m2,:] += overlap ./ ( im*pFreq.iwf .+ (                   - E1 + E2) )
								evalContributions[n1][4] += abs(overlap)
								evalContributions[n2][4] += abs(overlap)
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

	printGFnorm(gfdiagnorm, maximum( first.(evallist) )-E0)

	return gf_w, gf_iw, evalContributions
end 

####################################################################################################

"""
    getGF(evallist, eveclist, allstates, pModel, pSimulation, pFreq, pNumerics)

Evaluation of the Lehmann representation for the interacting finite-temperature Green's function
We sum over all Eigenstates `evallist` with electron number N and all other eigenstates n2 with electron number N-1
Since we have spin degeneracy, we only calculate the up GF, so we annihilate one up spin, so we also only sum over all S-1 states for given S(n1)
"""
function getGF( evallist::Array{Array{Eigenvalue,1},1},
               overlaps::Array{Complex{Float32},3},
					possibleTransitions::Array{Array{Array{Int64,1},1},1}, 
					NSperm::Array{Int64,1},
					pModel::ModelParameters,
					pSimulation::SimulationParameters,
					pFreq::FrequencyMeshes,
					pNumerics::NumericalParameters)::Tuple{SingleParticleFunction,SingleParticleFunction,Array{Array{Float64,1},1}}

	nflav=length(pSimulation.gf_flav)
	Nmax=UInt32(pModel.Nmax)
	beta = pSimulation.beta

	gfdiagnorm                    = zeros(Float64,nflav)                   # check normalization of GF
	gf_w::SingleParticleFunction  = zeros(nflav,nflav,length(pFreq.wf))
	gf_iw::SingleParticleFunction = zeros(nflav,nflav,length(pFreq.iwf))

	evalContributions::Array{Array{Float64,1},1} = []

	E0 = minimum( first.(evallist) )

	# prefill evalContribution array
	for n1=1:length(evallist)
		E1    = evallist[NSperm[n1]][1] -E0 
		N1    = round(Int64,evallist[NSperm[n1]][2])
		s1    = round(UInt64,evallist[NSperm[n1]][3])
		S1    = spinConfig(s1,N1,Nmax)
		push!( evalContributions, [N1,S1,E1,0.0] )
	end

	for n1=1:length(evallist)
	    if n1%(round(Int64,length(evallist)/10)) == 0
            print("\r"*lpad(round(Int64,n1*100.0/length(evallist)),4)*"%")
		end
		E1    = evallist[NSperm[n1]][1] -E0     # shift by E0 to avoid overflow
		N1    = round(Int64,evallist[NSperm[n1]][2])
		s1    = round(UInt64,evallist[NSperm[n1]][3])
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
				s2    = round(UInt64,evallist[NSperm[n2]][3])
				S2    = spinConfig(s2,N2,Nmax)

				expFac = exp(-beta*E1)+exp(-beta*E2)
	
				if ( expFac > pNumerics.cutoff )
					for b=a:length(pSimulation.gf_flav) # only upper triangular part
						if (abs(overlaps[2*b-1,n2,n1])>pNumerics.cutoff)

							ovrlp = overlaps[2*a-0,n1,n2]*overlaps[2*b-1,n2,n1]            # <n1|c_a|n2> * <n2|cdag_b|n1>
							gf_w[a,b,:]  += ovrlp ./ ( pFreq.wf     .+ (pNumerics.delta*im + E1 - E2) )
							gf_iw[a,b,:] += ovrlp ./ ( im*pFreq.iwf .+ (                   + E1 - E2) )
							evalContributions[n1][4] += abs(ovrlp)
							evalContributions[n2][4] += abs(ovrlp)
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
	for m1=1:nflav
		gf_iw[m1,m1,:] ./= gfdiagnorm[m1]
	end

	printGFnorm(gfdiagnorm, maximum( first.(evallist) )-E0)

	return gf_w, gf_iw, evalContributions
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
					pModel::ModelParameters,
					pSimulation::SimulationParameters,
					pFreq::FrequencyMeshes,
					pNumerics::NumericalParameters)::Tuple{TwoParticleFunction,Array{Array{Float64,1},1}}

	# we restrict this calculation to one orbital !
	#norb=1
	Nmax=UInt32(pModel.Nmax)
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

	evalContributions::Array{Array{Float64,1},1} = []
	# prefill the evalContributions array
	for m=1:length(evallist)
		Em    = evallist[NSperm[m]][1] -E0 
		Nm    = round(Int64,evallist[NSperm[m]][2])
		sm    = round(UInt64,evallist[NSperm[m]][3])
		Sm    = spinConfig(sm,Nm,Nmax)
		push!( evalContributions, [Nm,Sm,Em,0.0] )
	end

	# Now calculate the 2particle GF
	for m=1:length(evallist)
		if m%(max(1,round(Int64,length(evallist)/100.0))) == 0
            print("\r"*lpad(round(Int64,m*100.0/length(evallist)),4)*"%")
		end

		Em    = evallist[NSperm[m]][1] -E0     # shift by E0 to avoid overflow
		Nm    = round(Int64,evallist[NSperm[m]][2])
		sm    = round(UInt64,evallist[NSperm[m]][3])
		evecm = eveclist[NSperm[m]]
		Sm    = spinConfig(sm,Nm,Nmax)
		expFm = exp(-beta*Em) # exponential Boltzmann factor

		for n=1:length(evallist)
			En    = evallist[NSperm[n]][1] -E0
			Nn    = round(Int64,evallist[NSperm[n]][2])   
			sn    = round(UInt64,evallist[NSperm[n]][3])
			evecn = eveclist[NSperm[n]]
			Sn    = spinConfig(sn,Nn,Nmax)
			expFn = exp(-beta*En)

			for o=1:length(evallist)
				Eo    = evallist[NSperm[o]][1] -E0
				No    = round(Int64,evallist[NSperm[o]][2])   
				so    = round(UInt64,evallist[NSperm[o]][3])
				eveco = eveclist[NSperm[o]]
				So    = spinConfig(so,No,Nmax)
				expFo = exp(-beta*Eo)

				for p=1:length(evallist)
					Ep    = evallist[NSperm[p]][1] -E0
					Np    = round(Int64,evallist[NSperm[p]][2])   
					sp    = round(UInt64,evallist[NSperm[p]][3])
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
										evalContributions[m][4] += abs(overlap)
										evalContributions[n][4] += abs(overlap)
										evalContributions[o][4] += abs(overlap)
										evalContributions[p][4] += abs(overlap)
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
										evalContributions[m][4] += abs(overlap)
										evalContributions[n][4] += abs(overlap)
										evalContributions[o][4] += abs(overlap)
										evalContributions[p][4] += abs(overlap)
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
										evalContributions[m][4] += abs(overlap)
										evalContributions[n][4] += abs(overlap)
										evalContributions[o][4] += abs(overlap)
										evalContributions[p][4] += abs(overlap)
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
										evalContributions[m][4] += abs(overlap)
										evalContributions[n][4] += abs(overlap)
										evalContributions[o][4] += abs(overlap)
										evalContributions[p][4] += abs(overlap)
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
										evalContributions[m][4] += abs(overlap)
										evalContributions[n][4] += abs(overlap)
										evalContributions[o][4] += abs(overlap)
										evalContributions[p][4] += abs(overlap)
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
										evalContributions[m][4] += abs(overlap)
										evalContributions[n][4] += abs(overlap)
										evalContributions[o][4] += abs(overlap)
										evalContributions[p][4] += abs(overlap)
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

"""
    getPossibleTransitions(evallist,eveclist,allstates,flavors,pNumerics)

Calculate the overlap elements between all Eigenstates acting on c/c^dagger specified by flavors.
The function returns the overlap array and possibleTransitions for cdagger,c for 
the flavors specified in the flavors array (1=orb1,up, 2=orb1,dn, 3=orb2,up, ....)
"""
function getPossibleTransitions(evallist::Array{Array{Eigenvalue,1},1},
                                eveclist::Array{Eigenvector,1},
				   					  allstates::NSstates,
										  flavors::Array{Int64,1},
										  NSperm::Array{Int64,1},
										  pNumerics::NumericalParameters)
	Nmax = getNmaxFromAllstates(allstates)
	#NSperm = getNSperm(evallist)

	nstates = length(eveclist)
	nflavors = length(flavors)

	#upflavors = filter(x->isodd(x), flavors) # 1,3,5, ... are up spins
	#dnflavors = filter(x->iseven(x), flavors) # 2,4,6, ... are dn spins

	possibleTransitions = [ [ Int64[] for i=1:nstates] for a=1:nflavors*2 ] # TODO: Array like ['crea','up',m]?
	overlaps = zeros(Complex{Float32},nflavors*2,nstates,nstates)
	# we always calculate crea and anni overlap for each flavor, i.e. norb*2 
	# ordering is crea_1,anni_1, crea_2,anni2, ....

	# HOWTO define a CAmatrix without initializing?
	creamat = CAmatrix[ sparse([1],[1],[1]) for m in 1:nflavors]
	annimat = CAmatrix[ sparse([1],[1],[1]) for m in 1:nflavors]

	NScrea   = [ Int64[-1,-1,-1,-1] for m in 1:nflavors]
	NSanni   = [ Int64[-1,-1,-1,-1] for m in 1:nflavors]

	for n1=1:nstates
		if n1%(max(1,round(Int64,nstates/100.0))) == 0
            print("\r"*lpad(round(Int64,n1*100.0/nstates),4)*"%")
		end
		N1    = round(Int64,evallist[NSperm[n1]][2])
		s1    = round(UInt64,evallist[NSperm[n1]][3])
		S1    = spinConfig(s1,N1,Nmax)
		evec1 = eveclist[NSperm[n1]]

		for n2=1:nstates
			N2    = round(Int64,evallist[NSperm[n2]][2])
			s2    = round(UInt64,evallist[NSperm[n2]][3])
			S2    = spinConfig(s2,N2,Nmax)
			evec2 = eveclist[NSperm[n2]]

			for i=1:nflavors
				dS = 2*(flavors[i]%2)-1 # spin increases by 1 for crea if up flavor, etc

				# Cdagger  <1|c^dag|2>
				if N1==N2+1 && S1==S2+dS
					if NScrea[i]!=[N1,S1,N2,S2]    # do we need to recalculate the cmatrix?
						NScrea[i]=[N1,S1,N2,S2]
						creamat[i] = getCdagmatrix(flavors[i], allstates[N1+1][s1], allstates[N2+1][s2])
					end
					overlaps[2*i-1,n1,n2] = dot( evec1, creamat[i] * evec2 )
					if abs(overlaps[2*i-1,n1,n2]) > pNumerics.cutoff
						push!( possibleTransitions[2*i-1][n2], n1 )
					end
				end # if

				# C  <1|c|2>
				if N1==N2-1 && S1==S2-dS
					if NSanni[i]!=[N1,S1,N2,S2]    # do we need to recalculate the cmatrix?
						NSanni[i]=[N1,S1,N2,S2]
						annimat[i] = getCmatrix(flavors[i], allstates[N1+1][s1], allstates[N2+1][s2])
					end
					overlaps[2*i-0,n1,n2] = dot( evec1, annimat[i] * evec2 )
					if abs(overlaps[2*i-0,n1,n2]) > pNumerics.cutoff
						push!( possibleTransitions[2*i-0][n2], n1 )
					end
				end # if

			end # flavors

		end # n2
		#println( n1,": N=",N1,", S=",S1," -> ",length(possibleTransitions[1][n1])," (",round(length(possibleTransitions[1][n1])*100.0/nstates,digits=2)," %)" )
	end # n1
    println("\rdone!")
	return overlaps, possibleTransitions
end

#########################################################################################

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
