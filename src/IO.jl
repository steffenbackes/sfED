
"""
    writeMatrixGnuplot(filename, matrix)

Writes a Matrix `matrix(norb1,norb2)` in human readable (fixed column width format) form into a file.
Only writes the real part and writes a linebreak after each row so that Gnuplot sp w pm3d works
.
"""
function writeMatrixGnuplot(filename::String,matrix)
	
	norb1=size(matrix,1)
	norb2=size(matrix,2)

	open(filename, "w") do outf
		for m1=1:norb1
			for m2=1:norb2
				re = real(matrix[m1,m2])
#				im = imag(matrix[m1,m2])
				@printf(outf, "%i %i %.7E \n", m1,m2,re)
			end # m2
		write(outf, "\n" )
		end # m1
	end # close
end

"""
    writeGF(filename, G, w)

Writes a Green's function `G(norb,norb,nw)` in human readable (fixed column width format) form into a file.
`mgrid` is the associated Matsuabra grid (as indices or evaluated).
"""
function writeGF(filename::String,G::SingleParticleFunction,mgrid)
    @assert length(mgrid) <= size(G,3)
	
	norb=size(G,1)

	open(filename, "w") do outf
		for n=1:length(mgrid)
			ww = mgrid[n]
			@printf(outf, "%.7E   ", ww)
			for m1=1:norb
				for m2=1:norb
					reg = real(G[m1,m2,n])
					img = imag(G[m1,m2,n])
					@printf(outf, "%.7E  %.7E   ", reg, img)
				end # m2
			end # m1
			write(outf, "\n" )
		end # n
	end # close
end

"""
    writeGF2part(filename, G, w)

Writes a Green's function `G(nw,nw,nw)` in human readable (fixed column width format) form into a file.
`mgrid` is the associated Matsuabra grid (as indices or evaluated).
"""
function writeGF2part(filename::String,G::TwoParticleFunction,mgrid)
	nw1=size(G,1)
	nw2=size(G,2)
	nw3=size(G,3)

	open(filename, "w") do outf
		for n1=1:nw1
			ww1 = mgrid[n1]
			for n2=1:nw2
				ww1 = mgrid[n1]
				ww2 = mgrid[n2]
				@printf(outf, "%.7E  %.7E  ", ww1,ww2)

				reg = real(G[n1,n2,1])
				img = imag(G[n1,n2,1])
				@printf(outf, "%.7E  %.7E   ", reg, img)

				write(outf, "\n" )
			end # n2
			write(outf, "\n" )
		end # n1
	end # close
end
"""
    printStateInfo(allstates)

Print all states and N S quantum numbers
"""
function printStateInfo(allstates::NSstates)
	println("We have the following states:")
	Nmax = getNmaxFromAllstates(allstates)
	for n=0:Nmax
		for s=1:noSpinConfig(n,Nmax)
			println("N=",n,", S=",spinConfig(s,n,Nmax))
			for state in allstates[n+1][s]
				println(state)
			end
		end
	end
end

"""
    printEvalInfo(evallist, eveclist, allstates)

Print the Eigenvalues and some info
"""
function printEvalInfo(evallist::Array{Array{Eigenvalue,1},1},
                       eveclist::Array{Eigenvector,1},
   					     allstates::NSstates)
	Nmax = getNmaxFromAllstates(allstates)
	E0 = minimum(first.(evallist))
	perm = sortperm(first.(evallist))

	printlimit = round(Int64,min(10, length(evallist)/2))

	println("Eigenstates: (",printlimit," lowest and highest calculated)")
	for i=vcat(1:printlimit, length(evallist)-printlimit+1:length(evallist))
		E = evallist[perm][i][1] -E0
		N = round(Int64,evallist[perm][i][2])
		s = round(UInt64,evallist[perm][i][3])
		S = spinConfig(s,N,Nmax)
		evstate = eveclist[perm][i]
		permV = sortperm( evstate, by=abs, rev=true )
		@printf("E=%+10.5f, N=%3i, S=%+3i : ", E,N,S)
		# print the three largest contributions to the Eigenvector
		for j=1:min(2,length(evstate))
			@printf( "%4.2f", abs(evstate[permV[j]])^2 ) 
			print( "x",allstates[N+1][s][permV[j]],"  ")
		end
		println(" ")
		if (i==printlimit)
			println(".")
			println(".")
			println(".")
		end
	end
end



"""
    writeEvalContributionsSectors(filename, evalContributions)

write the weight contribtions of each Eigenstate sorted by N,S quantum numbers
"""
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
			E = evalContributions[perm[n1]][3]
			if writelabel == 1
				write(outf, "$n1  $val $E ($N,$S) \n" )
			else
				write(outf, "$n1  $val $E \n" )
			end
		end
	end # close
end

"""
    writeEvalContributions(file, evalContributions)

write the weight contribtions of each Eigenstate sorted by Eigenenergy
"""
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

"""
    printGFnorm(gfdiagnorm)

Print the calculated normalization of the diagonal Green's function elements
"""
function printGFnorm(gfdiagnorm::Array{Float64,1}, highestEval::Eigenvalue)
	println("Obtained the following normalization of the Green's function for all orbitals:")
	for n in gfdiagnorm
		@printf( "%5.3f ", n ) 
	end
	println(" ")
	if any(gfdiagnorm .< 0.95)
		@printf("Only %3i%% of all states are in the energy window [%+4.1f,%+4.1f] \n",minimum(gfdiagnorm)*100,-highestEval,highestEval)
		println("!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
		println("!!! Green's function is not normalized, increase Energy cutoff !!!")
		println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
	else
		@printf("%5.1f%% of all states are in the energy window [%+4.1f,%+4.1f] \n",minimum(gfdiagnorm)*100,-highestEval,highestEval)
	end
end

####################################################################################

"""
    writePossibleTransitions(filename,overlaps)

Write out the overlap elements between all Eigenstates acting on c/c^dagger
"""
function writePossibleTransitions(filename::String,overlaps::Array{Float32,2} )
	nstates = size(overlaps)[1]
	open(filename, "w") do outf
		for n1 in 1:nstates
			for n2 in 1:nstates
				@printf(outf, "%i  %i  %.7E \n", n1,n2,  overlaps[n1,n2])
				@printf(outf, "%i  %i  %.7E \n", n1,n2+1,overlaps[n1,n2])
			end # n2
			@printf(outf,"\n")
			for n2 in 1:nstates
				@printf(outf, "%i  %i  %.7E \n", n1+1,n2,  overlaps[n1,n2])
				@printf(outf, "%i  %i  %.7E \n", n1+1,n2+1,overlaps[n1,n2])
			end # n2
			@printf(outf,"\n")
		end # n1
	end # close file
end

################################################################################################################
