
# ============================ Comments ============================
# TODO: write macro that wraps all func(io::IO,...) in func(fname::String,...) and func(...) (stdou)
#
#
# ======================= Auxilliary Function =======================

# ========================= Main Functions =========================
"""
    writeMatrixGnuplot(filename, matrix)

Writes a Matrix `matrix(norb1,norb2)` in human readable (fixed column width format) form into a file.
Only writes the real part and writes a linebreak after each row so that Gnuplot sp w pm3d works
.
"""
function writeMatrixGnuplot(io::IO, matrix)
    norb1=size(matrix,1)
    norb2=size(matrix,2)
    for m1=1:norb1
        for m2=1:norb2
            re = real(matrix[m1,m2])
#           im = imag(matrix[m1,m2])
            @printf(io, "%i %i %.7E \n", m1,m2,re)
        end # m2
        write(io, "\n")
    end # m1
end

"""
    writeGF(filename, G, w)

Writes a Green's function `G(norb,norb,nw)` in human readable (fixed column width format) form into a file.
`mgrid` is the associated Matsuabra grid (as indices or evaluated).
"""
function writeGF(io::IO, G::SingleParticleFunction,mgrid)
    length(mgrid) <= size(G,3) || throw(DomainError(size(G,3) - length(mgrid), "size(G,3) must be larger than length of mgrid")) 
   
    norb=size(G,1)
    for n=1:length(mgrid)
        ww = mgrid[n]
        @printf(io, "%.7E   ", ww)
        for m1=1:norb
            for m2=1:norb
                reg = real(G[m1,m2,n])
                img = imag(G[m1,m2,n])
                @printf(io, "%.7E  %.7E   ", reg, img)
            end # m2
        end # m1
        write(io, "\n" )
    end # n
end

"""
    writeGF2part(filename, G, w)

Writes a Green's function `G(nw,nw,nw)` in human readable (fixed column width format) form into a file.
`mgrid` is the associated Matsuabra grid (as indices or evaluated).
"""
function writeGF2part(io::IO, G::TwoParticleFunction,mgrid)
    nw1=size(G,1)
    nw2=size(G,2)
    nw3=size(G,3)

    for n1=1:nw1
        ww1 = mgrid[n1]
        for n2=1:nw2
            ww1 = mgrid[n1]
            ww2 = mgrid[n2]
            @printf(io, "%.7E  %.7E  ", ww1,ww2)

            reg = real(G[n1,n2,1])
            img = imag(G[n1,n2,1])
            @printf(io, "%.7E  %.7E   ", reg, img)

            write(io, "\n" )
         end # n2
         write(io, "\n" )
    end # n1
end
"""
    writeStateInfo(allstates)

Write all states and N S quantum numbers
"""
function writeStateInfo(io::IO, allstates::NSstates)
   println(io, "We have the following states:")
   Nmax = getNmaxFromAllstates(allstates)
   for n=0:Nmax
      for s=1:noSpinConfig(n,Nmax)
         println(io, "N=",n,", S=",spinConfig(s,n,Nmax))
         for state in allstates[n+1][s]
            println(io, state)
         end
      end
   end
end

"""
    writeEvalInfo(evallist, eveclist, allstates)

Write the Eigenvalues and some info
"""
function writeEvalInfo(io::IO, eigenspace::Eigenspace,
                               fockstates::Fockstates)
   Nmax = fockstates.Nmax
   Nstates = eigenspace.Nstates
   E0 = minimum(eigenspace.evals)
   perm = sortperm(eigenspace.evals) # sort by Eigenenergy

   printlimit = round(Int64,min(10, Nstates/2))

   println(io, "Eigenstates: (",printlimit," lowest and highest)")
   for i=vcat(1:printlimit, Nstates-printlimit+1:Nstates)
      E,evstate = eigenspace[perm[i]]
      E -= E0
      N = eigenspace.Nel[perm[i]]
      S = eigenspace.Spin[perm[i]]
      is = indexSpinConfig(S,N,Nmax)

      permV = sortperm( evstate, by=abs, rev=true )
      @printf(io, "E=%+10.5f, N=%3i, S=%+3i : ", E,N,S)
      # print the three largest contributions to the Eigenvector
      for j=1:min(2,length(evstate))
         @printf(io,  "%4.2f", abs(evstate[permV[j]])^2 ) 
         print(io, "x",fockstates[N,is,permV[j]],"  ")
      end
      println(io, " ")
      (i==printlimit) && println(io, ".\n.\n.")
   end
end



"""
    writeEvalContributionsSectors(filename, evalContributions)

write the weight contribtions of each Eigenstate sorted by N,S quantum numbers
"""
function writeEvalContributionsSectors(io::IO, evalContributions::Array{Float64,2})
   #perm = sortperm([ evalContributions[:,2:end][i,:] for i in 1:size(evalContributions[:,2:end])[1]] )                   # sort by N,S, then Eval
   NSlist = []
   for i=1:size(evalContributions)[1]
      push!( NSlist, evalContributions[i,:])
   end
   perm = sortperm(NSlist)

    for n1=1:size(evalContributions)[1]
        writelabel = 1
        if n1>1
            if round.(Int64, evalContributions[perm[n1-1],1:2]) != round.(Int64, evalContributions[perm[n1],1:2] )
                write(io, "\n" )
                writelabel = 1
            else
                writelabel = 0
            end
        end
   
        val = evalContributions[perm[n1],4]
        N = round(Int64, evalContributions[perm[n1],1])
        S = round(Int64, evalContributions[perm[n1],2])
        E = evalContributions[perm[n1],3]
        if writelabel == 1
            write(io, "$n1  $val $E ($N,$S) \n" )
        else
            write(io, "$n1  $val $E \n" )
        end
    end
end

"""
    writeEvalContributions(file, evalContributions)

write the weight contribtions of each Eigenstate sorted by Eigenenergy
"""
function writeEvalContributions(io::IO, evalContributions::Array{Float64,2})
    # sort the contributions according to energy
    evals = []
    for n1=1:size(evalContributions)[1]
        push!(evals, evalContributions[n1,3] )
    end
    perm = sortperm(evals)
    for n1=1:size(evalContributions)[1]
        val = evalContributions[perm[n1],4]
        E = evalContributions[perm[n1],3]
        write(io, "$n1 $E  $val  \n" )
    end
end

"""
    writeGFnorm(gfdiagnorm)

Write the calculated normalization of the diagonal Green's function elements
"""
function writeGFnorm(io::IO, gfdiagnorm::Array{Float64,1})
   println("Obtained the following normalization of the Green's function for all orbitals:")
   for n in gfdiagnorm
      @printf(io,  "%5.3f ", n )
   end
   println(" ")
   if any(gfdiagnorm .< 0.95)
      #@printf(io, "Only %3i%% of all states are in the energy window [%+4.1f,%+4.1f] \n",minimum(gfdiagnorm)*100,-highestEval,highestEval)
      @printf(io, "Only %3i%% of all states obtained! \n",minimum(gfdiagnorm)*100)
      println("!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      println("!!! Green's function is not normalized !!!")
      println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
   else
      @printf(io, "%5.1f%% of all states obtained. \n",minimum(gfdiagnorm)*100)
   end
end

"""
    writeTransitionsOverlaps(filename,overlaps)

Write out the overlap elements between all Eigenstates acting on c/c^dagger
"""
function writeTransitionsOverlaps(io::IO, overlaps::Array{Complex{Float32},3} )
   nstates = size(overlaps)[2]
   for n1 in 1:nstates
     for n2 in 1:nstates
        @printf(outf, "%i  %i  ", n1,n2)
        for val in overlaps[:,n1,n2]
           @printf(outf, "%.7E ", abs(val))
        end
        @printf(outf, "\n")

        @printf(outf, "%i  %i  ", n1,n2+1)
        for val in overlaps[:,n1,n2]
           @printf(outf, "%.7E ", abs(val))
        end
        @printf(outf, "\n")
     end # n2
     @printf(outf,"\n")
     for n2 in 1:nstates
        @printf(outf, "%i  %i  ", n1+1,n2)
        for val in overlaps[:,n1,n2]
           @printf(outf, "%.7E ", abs(val))
        end
        @printf(outf, "\n")

        @printf(outf, "%i  %i  ", n1+1,n2+1)
        for val in overlaps[:,n1,n2]
           @printf(outf, "%.7E ", abs(val))
        end
        @printf(outf, "\n")
     end # n2
     @printf(outf,"\n")
     end # n1
end


# ========================= Boilerplate Generator =========================
io_functions = (:writeMatrixGnuplot, :writeGF, :writeGF2part, :writeStateInfo, :writeEvalInfo, :writeEvalContributionsSectors, :writeEvalContributions, :writeGFnorm, :writeTransitionsOverlaps)

# For all write functions f(io::IO, args...) generate wrappers f(args...) which writes
# directly to stdout and f(s::String, args...) which opens/closes a file and writes to it.
for ff in io_functions
    @eval $ff(s::String, arg, args...) = (open(s,"w") do io; $ff(io,arg,args...); end)
    @eval $ff(arg, args...) = $ff(stdout, arg, args...)
end
