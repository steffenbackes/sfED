
########################################################################
# Writes a Green's function G(norb,norb,nw) like object into a file
function writeGF(filename::String,
                       gf::Array{Complex{Float64},3},
                        w)            # w can be linspace or Float64 array
	open(filename, "w") do outf
		for n=1:length(w)
			ww = w[n]
			write(outf, "$ww  " )
			for m1=1:norb
				for m2=1:norb
					reg = real(gf[m1,m2,n])
					img = imag(gf[m1,m2,n])
					write(outf, "$reg  $img  " )  # does one really have to do it like this??
				end # m2
			end # m1
			write(outf, "\n" )
		end # n
	end # close
end
########################################################################

#######################################################################
# Print all states and N S quantum numbers #########################
function printStateInfo(allstates::Array{Array{Array{Array{Int64,1},1},1},1})
	println("We have the following states:")
	for n=0:Nmax
		for s=1:noSpinConfig(n,Nmax)
			println("N=",n,", S=",spinConfig(s,n,Nmax))
			for state in allstates[n+1][s]
				println(state)
			end
		end
	end
end
#######################################################################

#######################################################################
# Print the Eigenvalues and some info
function printEvalInfo( evallist::Array{Array{Float64,1},1},
                        eveclist::Array{Array{Complex{Float64},1},1},
							  allstates::Array{Array{Array{Array{Int64,1},1},1},1})
	println("Eigenstates:")
	for i=1:length(evallist)
		E = evallist[i][1]
		N = round(Int64,evallist[i][2])
		s = round(Int64,evallist[i][3])
		S = spinConfig(s,N,Nmax)
		evstate = eveclist[i]
		perm = sortperm( evstate, by=abs, rev=true )
		@printf("E=%+10.5f, N=%3i, S=%+3i : ", E,N,S)
		# print the three largest contributions to the Eigenvector
		for j=1:min(2,length(evstate))
			@printf( "%4.2f", abs(evstate[perm[j]])^2 ) 
			print( "x",allstates[N+1][s][perm[j]],"  ")
		end
		println(" ")
	end
end
#######################################################################


#######################################################################
# write the weight contribtions of each Eigenstate sorted by N,S quantum numbers
function writeEvalContributionsSectors(filename::String, evalContributions::Array{Array{Float64,1},1})
	perm = sortperm(evalContributions)                   # sort by N,S, then Eval
	open(filename, "w") do outf
		for n1=1:length(evalContributions)
			writelabel = 1
			if n1>1
				if round.(Int64, evalContributions[perm[n1-1]][1:2]) != round.(Int64, evalContributions[perm[n1]][1:2] )
					write(outf, "\n" )
					writelabel = 1
				else
					writelabel = 0
				end
			end
	
			val = evalContributions[perm[n1]][4]
			N = round(Int64, evalContributions[perm[n1]][1])
			S = round(Int64, evalContributions[perm[n1]][2])
			if writelabel ==1
				write(outf, "$n1  $val ($N,$S) \n" )
			else
				write(outf, "$n1  $val  \n" )
			end
		end
	end # close
end
#######################################################################

#######################################################################
# write the weight contribtions of each Eigenstate sorted by Eigenenergy
function writeEvalContributions(filename::String, evalContributions::Array{Array{Float64,1},1})
	# sort the contributions according to energy
	evals = []
	for n1=1:length(evalContributions)
		push!(evals, evalContributions[n1][3] )
	end
	perm = sortperm(evals)
	open(filename, "w") do outf
		for n1=1:length(evalContributions)
			val = evalContributions[perm[n1]][4]
			E = evalContributions[perm[n1]][3]
			write(outf, "$n1 $E  $val  \n" )
		end
	end # close
end
#######################################################################
