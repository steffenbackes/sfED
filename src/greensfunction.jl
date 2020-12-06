"""
    getG0(eps, tmatrix, w)

Construct the noninteracting Green's function from onsite energy `eps`,
hoppings between sites `i` and `j`, `tmatrix` and Matsubara grid `w`.
"""
function getG0(eps::Array{Float64,1},tmatrix::Array{Float64,2},
               
					      w )
	gf0 = zeros(ComplexF64,norb,norb,nw)
	for n=1:nw
		gf0[:,:,n] = inv( I*(w[n] + mu) - Diagonal(eps[1:2:end]) - tmatrix ) # eps is defined for spins, so only take orbital part
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
                    subspace2::Array{Array{Int64,1},1}, 
						  subspace1::Array{Array{Int64,1},1})
	dim1 = length(subspace1)
	dim2 = length(subspace2)

	cmatrix = zeros(Int64,dim2,dim1)
	# now annihilate the particle in all basis states and find the corresponding state in the N-1,S-1 space
	for j=1:dim1
		if subspace1[j][anni] == 1   # if there is a particle to annihilate..
			state = copy( subspace1[j] )
			state[anni] = 0
			a1sgn = (-1)^sum(state[1:anni]) # count the particles before anni

			# now find this state in the smaller subspace
			i = findfirst(map(x-> all(x .== state), subspace2))
			cmatrix[i,j] = a1sgn
		end
	end # j
	return cmatrix
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
                       subspace1::Array{Array{Int64,1},1}, 
							  subspace2::Array{Array{Int64,1},1})
	dim1 = length(subspace1)
	dim2 = length(subspace2)

	cdagmatrix = zeros(Int64,dim1,dim2)
	# now create the particle in all basis states and find the corresponding state in the N,S space
	for j=1:dim2
		if subspace2[j][crea] == 0   # if there is space to create one particle...
			state = copy( subspace2[j] )
			c1sgn = (-1)^sum(state[1:crea]) # count the particles before crea
			state[crea] = 1

			# now find this state in the N,S subspace
			i = findfirst(map(x-> all(x .== state), subspace1))
			cdagmatrix[i,j] = c1sgn
		end
	end # j
	return cdagmatrix
end
#######################################################################

"""
    getGF(evallist, eveclist, allstates)

Evaluation of the Lehmann representation for the interacting finite-temperature Green's function
We sum over all Eigenstates `evallist` with electron number N and all other eigenstates n2 with electron number N-1
Since we have spin degeneracy, we only calculate the up GF, so we annihilate one up spin, so we also only sum over all S-1 states for given S(n1)
"""
function getGF( evallist::Array{Array{Float64,1},1},
                eveclist::Array{Array{Complex{Float64},1},1},
					allstates::Array{Array{Array{Array{Int64,1},1},1},1})
	NSvalues   = Int64[-1,-1,-1,-1]
	cmatrix    = zeros(Int64,norb,1,1)
	cdagmatrix = zeros(Int64,norb,1,1)
	gf_w       = zeros(ComplexF64,norb,norb,nw)
	gf_iw      = zeros(ComplexF64,norb,norb,nw)
	evalContributions::Array{Array{Float64,1},1} = []

	# First create a sortting for all Eigenstates in order of increasing N,S, which makes the GF generation more efficient
	NSperm = getNSperm(evallist)

	E0 = minimum( first.(evallist) )

	for n1=1:length(evallist)
	    if n1%(Int64(length(evallist)/10)) == 0
            print("\r"*lpad(Int64(n1*100.0/length(evallist)),4)*"%")
		end

		E1    = evallist[NSperm[n1]][1] -E0     # shift by E0 to avoid overflow
		N1    = round(Int64,evallist[NSperm[n1]][2])
		s1    = round(Int64,evallist[NSperm[n1]][3])
		evec1 = eveclist[NSperm[n1]]
		S1    = spinConfig(s1,N1,Nmax)
		push!( evalContributions, [N1,S1,E1,0.0] )
	
		for n2=1:length(evallist)
			E2    = evallist[NSperm[n2]][1] -E0
			N2    = round(Int64,evallist[NSperm[n2]][2])   
			s2    = round(Int64,evallist[NSperm[n2]][3])
			evec2 = eveclist[NSperm[n2]]
			S2    = spinConfig(s2,N2,Nmax)
	
			# Exclude transitions too high in energy or all Eigenstates 2 which are not N-1,S-1 
			if ( (exp(-beta*E1)+exp(-beta*E2))>cutoff && N2==N1-1 && S2==S1-1 )
	
				# If we have not been dealing with this N,S combination in the loop before, we need to generate the right Cmatrix and Cdagmatrix
				if ( NSvalues!=[N1,N2,S1,S2] )
					NSvalues = [N1,N2,S1,S2]
	
					# we need to generate cmatrix and cdagmatrix for all orbitals
					dim1 = length(allstates[N1+1][s1])
					dim2 = length(allstates[N2+1][s2])
					cmatrix    = zeros(Int64,norb,dim2,dim1)  # go from subspace 1 to 2 by annihilation
					cdagmatrix = zeros(Int64,norb,dim1,dim2)  # go from subspace 2 to 1 by creation
					for m1=1:norb
						c1 = 2*m1-1
						a1 = 2*m1-1
						cmatrix[m1,:,:]    =    getCmatrix(a1, allstates[N2+1][s2], allstates[N1+1][s1])
						cdagmatrix[m1,:,:] = getCdagmatrix(c1, allstates[N1+1][s1], allstates[N2+1][s2])
					end # m1 loop
				end # if we need to update cmatrix,cdagmatrix
	
				#Now we have the proper N,N-1 state and the c and c^dagger matrices, so evaluate the Lehmann representation
				# Obtain all matrix elements
				for m1=1:norb
					# first the diagonal elements
					overlap =  dot( evec1, cdagmatrix[m1,:,:]*evec2 ) * dot( evec2, cmatrix[m1,:,:]*evec1 )
					if abs(overlap)>cutoff
						gf_w[m1,m1,:]  += overlap * (exp(-beta*E1)+exp(-beta*E2)) ./ ( LinRange(wmin,wmax,nw) .+ (delta*im - E1 + E2) )
						gf_iw[m1,m1,:] += overlap * (exp(-beta*E1)+exp(-beta*E2)) ./ ( ( 2 .* collect(1:nw) .-1).*(im*pi/beta) .+ (- E1 + E2) )
						evalContributions[n1][4] += abs(overlap * (exp(-beta*E1)+exp(-beta*E2)))
						evalContributions[n2][4] += abs(overlap * (exp(-beta*E1)+exp(-beta*E2)))
					end
	
					# then upper offdiagonal
					for m2=m1+1:norb
						overlap =  dot( evec1, cdagmatrix[m1,:,:]*evec2 ) * dot( evec2, cmatrix[m2,:,:]*evec1 )
						if abs(overlap)>cutoff
							gf_w[m1,m2,:]  += overlap * (exp(-beta*E1)+exp(-beta*E2)) ./ ( LinRange(wmin,wmax,nw) .+ (delta*im - E1 + E2) )
							gf_iw[m1,m2,:] += overlap * (exp(-beta*E1)+exp(-beta*E2)) ./ ( ( 2 .* collect(1:nw) .-1).*(im*pi/beta) .+ (- E1 + E2) )
							evalContributions[n1][4] += abs(overlap * (exp(-beta*E1)+exp(-beta*E2)))
							evalContributions[n2][4] += abs(overlap * (exp(-beta*E1)+exp(-beta*E2)))
						end
					end # m2 loop
				end # m1 loop
	
			end # if exp(-beta*E) > cutoff
		end # n2 loop
	end # n1 loop
    println("\rdone!")
	
	# copy the offdiagonals
	for m1=1:norb
		for m2=m1+1:norb
			gf_w[m2,m1,:] = gf_w[m1,m2,:]
			gf_iw[m2,m1,:] = gf_iw[m1,m2,:]
		end
	end
	
	#normalize
	gf_w ./= getZ(evallist)
	gf_iw ./= getZ(evallist)

	return gf_w, gf_iw, evalContributions
end 

"""
    getSigma(G0, G)

Calculate the Selfenergy from `G0` and `G`
"""
function getSigma(G0::Array{Complex{Float64},3}, G::Array{Complex{Float64},3})
	sigma = zeros(ComplexF64,size(G))
	for i=1:nw
		sigma[:,:,i]  = inv( G0[:,:,i]) - inv( G[:,:,i])
	end
	return sigma
end
