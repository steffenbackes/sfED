
########################################################################
# Return numnber of configurations with different total spin S possible for given n and nMax
function noSpinConfig(n::Int64, Nmax::Int64)::Int64
	round(Int64,  Nmax/2-abs(n-Nmax/2)+1 )
end 
########################################################################

########################################################################
# Return the i-th possible total spin value S for given n,nMax. i starts at 1 (eg S=-2,0,2) 
# Our electrons have spin +- 1 !!!!!!!!!!!!!!!!!!!
function spinConfig(i::Int64,n::Int64,Nmax::Int64)::Int64
	 round(Int64, -(Nmax/2-abs(n-Nmax/2)) + (i-1)*2 )
end
########################################################################

########################################################################
# Return the index i corresponding to a spin config S for given n,nMax
function indexSpinConfig(S::Int64,n::Int,Nmax::Int)::Int64
	round(Int64,  ( S + (Nmax/2-abs(n-Nmax/2)) )/2 +1 )
end
########################################################################


########################################################################
# this function returns the total spin S of a given state
function getSpin(state::Array{Int64,1})
	# Array to calculate the total spin [1,-1,1,-1, etc]
	spinsign = Int64.( ones(Nmax) )
	for i=2:2:Nmax
		spinsign[i] = Int64(-1)
	end
	return sum(spinsign .* state)
end
########################################################################

########################################################################
# Here we create the array that stores all states as an integer array sorted by N,S quantum numbers
function generateStates()
	# Create an empty list that will store all states, sorted for given particle number and spin
	allstates::Array{Array{Array{Array{Int64,1},1},1},1} = []
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
		state = parse.(Int64, split( bitstring(i-1)[end-Nmax+1:end],"") )
		N = sum(state)
		S = getSpin(state)
		iS = indexSpinConfig(S,N,Nmax)
		push!( allstates[N+1][iS], state )
	end
	return allstates
end
########################################################################



#######################################################################
# Return permutation that sorts the eigenvalue/vector array according to N,S quantum numbers (increasing order)
function getNSperm(evallist::Array{Array{Float64,1},1})
	NS=[]
	for i=1:length(evallist)
		push!(NS, [evallist[i][2],evallist[i][3]] )   # extract the n,s quantum number into an array
	end
	return sortperm(NS)
end
#######################################################################

