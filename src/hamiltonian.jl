const Hamiltonian = Matrix{ComplexF64} #SparseMatrixCSC{ComplexF64,Int64}   # matrix type for Hamiltonian matrix

include("MatrixTemplates.jl")

"""
    getEps(pNumerics,pModel)

Return the local orbital levels, define them by hand here. Values are perturbed by `smallValue`
TODO: this should be overloaded to accept text file or command line inputs.
"""
function getEps(fockstates::Fockstates, pNumerics::NumericalParameters)
   eps = zeros(Float64,fockstates.norb*2)
   for i=0:fockstates.norb-1              # Just shift one orbital down, the other up by +-1
      eps[2*i+1] = 0.0*(-1)^i    # up spin
      eps[2*i+2] = eps[2*i+1]    #dn spin
   end

   if fockstates.norb==5  # Use Julian's AIM benchmark system
      for s=1:2
      #beta=11
         eps[2*0+s] = 0.0
         eps[2*1+s] = 1.03583054366269
         eps[2*2+s] = 0.178882722915798  
         eps[2*3+s] = -1.03583054366269 
         eps[2*4+s] = -0.178882722915798 
         #beta = 10
         #eps[2*1+s] = 1.03530343377897
         #eps[2*2+s] = 0.175195925728476  
         #eps[2*3+s] = -1.03530343377897  
         #eps[2*4+s] = -0.175195925728476 
      end
   end

   # add a very small random term to each local level to lift degeneracy and improve numerical stability
#   for i=1:length(eps)
#      eps[i] += rand([-1,1]) * rand(Float64) * pNumerics.cutoff
#   end
   return eps
end


"""
    getHopping(istate, jstate, tmatrix)

return the hopping contribution for states `<i| |j>`
"""
function getHopping( istate::Fockstate,
                     jstate::Fockstate,
                     tmatrix::Matrix,
                     pNumerics::NumericalParameters)
    # we act c^dagger c on the state |j> and check overlap with <i|
    norb = size(tmatrix)[1]
    htmp = 0.0
    for m1=1:norb
        for m2=1:norb
            # no hopping for same orbitals (this is in eps)
            if ( abs(tmatrix[m1,m2]) > pNumerics.cutoff && m1!=m2 )
                tval = tmatrix[m1,m2]
                for s=0:1 # spin
                    a1 = (m1-1)*2+s  +1
                    c1 = (m2-1)*2+s  +1
                    htmp += tval*getMatrixElem1Particle(c1,a1,istate,jstate)
                end
            end # if t>cutoff
        end # m2
    end # m1
    return htmp
end


"""
    getUdensity(state, Umatrx, Jmatrix)

return the  density density part of the Coulomb interaction contribution
"""
function getUdensity(  state::Fockstate,
                     Umatrix::Matrix,
                     Jmatrix::Matrix,
                     pNumerics::NumericalParameters)
    norb = size(Umatrix)[1]
    htmp = 0.0
    for m1=1:norb
        n1up = state[2*m1-1]
        n1dn = state[2*m1-0]
        htmp += Umatrix[m1,m1]*n1up*n1dn
        for m2=m1+1:norb
            if abs(Umatrix[m1,m2]) > pNumerics.cutoff
                n2up = state[2*m2-1]
                n2dn = state[2*m2-0]
                htmp += Umatrix[m1,m2]*n1up*n2dn
                htmp += Umatrix[m1,m2]*n1dn*n2up
                htmp += (Umatrix[m1,m2]-Jmatrix[m1,m2])*n1up*n2up
                htmp += (Umatrix[m1,m2]-Jmatrix[m1,m2])*n1dn*n2dn
            end
        end # m2
        end # m1
    return htmp
end

"""
    getUnondensity(istate, jstate, Jmatrix)

return the spin-flip and pair-hopping part of the interaction contribution
"""
function getUnondensity( istate::Fockstate,
                         jstate::Fockstate,
                        Jmatrix::Matrix,
                        pNumerics::NumericalParameters)
    norb = size(Jmatrix)[1]
    htmp = 0.0
    for m1=1:norb
        for m2=1:norb
        if ( abs(Jmatrix[m1,m2]) > pNumerics.cutoff && m1!=m2 )
            Jval = Jmatrix[m1,m2]

            # spin flip !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            a4 = (m2-1)*2+0  +1 # +1 since Julia starts indexing at 1
            a3 = (m1-1)*2+1  +1
            c2 = (m2-1)*2+1  +1
            c1 = (m1-1)*2+0  +1

            htmp += Jval * getMatrixElem2Particle(c1,c2,a3,a4,istate,jstate)

            # pair hopping !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            a4 = (m2-1)*2+0  +1 # +1 since Julia starts indexing at 1
            a3 = (m2-1)*2+1  +1
            c2 = (m1-1)*2+1  +1
            c1 = (m1-1)*2+0  +1

            htmp += Jval * getMatrixElem2Particle(c1,c2,a3,a4,istate,jstate)

            end # if Jmat>cutoff
        end # m2
    end # m1
    return htmp
end

"""
    getMatrixElem1Particle(c1, a1, istate, jstate)

Calculate the matrix element for a one particle operator given the creation and annihilation operator
with indices `c1` and `a1` repsectively.
"""
function getMatrixElem1Particle(c1::Int64,
                                a1::Int64,
                         istate::Fockstate,
                         jstate::Fockstate)
    tmp = 0.0
    # check that there is a particle to annihilate and an empty state to create when going j->i,
    # and that the reverse holds for the <i| state
    # then we can compare the states by counting the difference between them, it has to be exactly
    # 2, since they differ only in a1 and c1 position. Then |i> will be <j| after applying c^dagger c
    if ( jstate[a1] * (1-jstate[c1]) * (1-istate[a1]) * istate[c1] ==1 && sum(istate .!= jstate)==2 )
        a1sgn = (-1)^sum(jstate[1:a1-1])  # We need to subtract 1 because there is a particle at a1 in |i>
        c1sgn = (-1)^sum(istate[1:c1-1])    # We need to subtract 1 because there is a particle at c1 in <j|
        tmp = a1sgn*c1sgn
    end # if

    return tmp
end


"""
    getMatrixElem2Particle(c1,c2,a3,a4,istate,jstate)

Calculate the matrix element for a two particle operator given the creation and annihilation operators
with indices `c1`, `c2` and `a3`,`a4`.
"""
function getMatrixElem2Particle(c1::Int64,c2::Int64,
                                a3::Int64,a4::Int64,
                                istate::Fockstate,
                                jstate::Fockstate)
    tmp = 0.0 
    if (jstate[a4]==1 && jstate[c2]==0  &&      # check that there are two particles to annihilate in |j>
        jstate[a3]==1 && jstate[c1]==0  &&      # and space for two to create in <i|
        istate[a4]==0 && istate[c2]==1  &&
        istate[a3]==0 && istate[c1]==1 )

        state = copy(jstate)
        state[a4] = 0                   # First destroy particle
        a4sgn = (-1)^sum(state[1:a4])   # Then we can safely count the #particles before the particle in the state
        state[a3] = 0                   # Destroy next particle
        a3sgn = (-1)^sum(state[1:a3])   # safely count the #particles
        c2sgn = (-1)^sum(state[1:c2])   # We can safely count the #particles before creating since state[c2]==0
        state[c2] = 1                   # Destroy next particle
        c1sgn = (-1)^sum(state[1:c1])   # We can safely count the #particles before creating since state[c2]==0
        state[c1] = 1                   # Destroy next particle

        tmp = (state==istate) ? a4sgn*a3sgn*c2sgn*c1sgn : 0.0
    end # if
    return tmp
end

"""
    getHamiltonian(eps,tmatrix,Umatrix,Jmatrix,states)

Constructs hamiltonian from onsite energies `eps`, hoppings `tmatrix`, coulomb interactions `Umatrix` and `Jmatrix`
over precomputed `states`.

TODO: moved from sparse to StaticArray. Write a decision routine and chose between static (small), normal, sparse (large) version for very large basis (norb > 12?)
"""
function getHamiltonian(eps::Vector{Float64},tmatrix::Matrix,
                        Umatrix::Matrix,Jmatrix::Matrix, 
                        mu::Float64,aim::Int64,
                        states::Vector{Fockstate},
                        pNumerics::NumericalParameters)::Hamiltonian

    Hsize = length(states)
    H = Matrix{ComplexF64}(undef, Hsize, Hsize)
    for i=1:Hsize
        # set the diagonals (chemical potential differs between AIM and Hubbard model)
        H[i,i]  = (aim == 1) ? -mu*sum(states[i][1:2]) : -mu*sum(states[i])
        H[i,i] += sum(eps .* states[i])           # onsite levels
        # get Density-Density interaction terms
        H[i,i] += getUdensity(states[i], Umatrix, Jmatrix, pNumerics)
        # now offdiagonals
        for j=i+1:Hsize
            # Hopping terms ###########################################
            Hijtmp  = getHopping(states[i], states[j], tmatrix,pNumerics)
            # Pair-hopping and spin-flip terms of the Coulomb interaction
            Hijtmp += getUnondensity(states[i], states[j], Jmatrix,pNumerics)

            H[i,j] = Hijtmp 
            H[j,i] = conj(Hijtmp)
        end # j ket
    end # i bra

    # Construct the sparse matrix
    return H
end


"""
    getEvalEvecs(Hamiltonian)

Diagonalize a given `Hamiltonian`. Eigenvalues will be cast to real, since the Hamiltonian is hermitian.
"""
function getEvalEvecs(hamiltonian::Hamiltonian)::Tuple{Vector{Eigenvalue}, EigenvectorMatrix}
      HamDense = Hermitian(Matrix(hamiltonian))
      evals,evecs = eigen(HamDense)
      return real(evals), evecs
end

"""
    getEvalveclist(eps,tmatrix,Umatrix,Jmatrix,mu,fockstates,pNumerics)

Create the N,S submatrices of the Hamiltonian, solve it and return the Eigenvalues and Vectors in a List
"""
function getEvalveclist(eps::Vector{Float64},tmatrix::Matrix,
                        Umatrix::Matrix,Jmatrix::Matrix,
                        mu::Float64, aim::Int64,
                        fockstates::Fockstates,
                        pNumerics::NumericalParameters)
    Nmax = fockstates.norb*2
    evallist = [[] for n=0:Nmax]
    eveclist = [[] for n=0:Nmax]

    for n=0:Nmax
        for s=1:noSpinConfig(n,Nmax)
            dim = length(fockstates.states[n+1][s])
            print("Constructing Hamiltonian(",dim,"x",dim,"), N=",n,", S=",spinConfig(s,n,Nmax),"... ")

            # now get the Hamiltonian submatrix spanned by all states <i| |j> in the N,S space (sparse matrix)
            hamiltonian = getHamiltonian(eps,tmatrix,Umatrix,Jmatrix,mu,aim,fockstates.states[n+1][s],pNumerics)
            els = prod(size(hamiltonian))

            # (n==8 && spinConfig(s,n,Nmax)==0) && (writeMatrixGnuplot("hamil.dat",Matrix(hamiltonian));  exit())

            print("Diagonalizing Hamiltonian... ")
            evals,evecs = getEvalEvecs(hamiltonian)

            push!(evallist[n+1], evals)
            push!(eveclist[n+1], [evecs[:,i] for i=1:length(evals)])
            println("Done!")
        end # s
    end # n
    return evallist,eveclist
end


"""
    getZ(evallist)

Calculate partition function from the eigenvalues in `evallist`.
"""
function getZ(eigenspace::Eigenspace, beta::Float64)
    evals = [ e for nse in eigenspace.evals for se in nse for e in se  ]
    return sum( exp.( -beta .*( evals .- eigenspace.E0) ) )
end

"""
    getN(eigenspace)

Calculate expectation value of occupation
"""
function getN(eigenspace::Eigenspace, beta::Float64, fockstates::Fockstates)
    Nmax = fockstates.norb*2
    Z = getZ(eigenspace, beta)
    n_ms = zeros(Float64,fockstates.norb*2)  # Filling per orbital and spin
    for n=0:Nmax
        for s=1:noSpinConfig(n,Nmax)
            dim = length(eigenspace.evals[n+1][s])
            for i=1:dim
                for j=1:dim
                    n_ms += abs(eigenspace.evecs[n+1][s][i][j])^2 .* fockstates.states[n+1][s][j] * exp(-beta*(eigenspace.evals[n+1][s][i]-eigenspace.E0))
                end # j
            end # i
        end # s
    end # n
    return n_ms/Z
end

"""
    getE(eigenspace)

Calculate expectation value of Energy
"""
function getE(eigenspace::Eigenspace, beta::Float64)
    evals = [ e for nse in eigenspace.evals for se in nse for e in se  ]
    return sum( evals .* exp.( -beta .*( evals .- eigenspace.E0) ) ) / getZ(eigenspace,beta)
end

"""
    getNN(eigenspace)

Calculate expectation value of all double occupation combinations
"""
function getNN(eigenspace::Eigenspace, beta::Float64, fockstates::Fockstates)::Array{Float64,3}
    Nmax = fockstates.norb*2
    Z = getZ(eigenspace, beta)
    nn = zeros(Float64,3,fockstates.norb,fockstates.norb)  # Nup*Ndn, Nup*Nup, Ndn*Ndn as orbital matrix
    for n=0:Nmax
        for s=1:noSpinConfig(n,Nmax)
            dim = length(eigenspace.evals[n+1][s])
            for i=1:dim
                for j=1:dim
                    weight = abs(eigenspace.evecs[n+1][s][i][j])^2 * exp(-beta*(eigenspace.evals[n+1][s][i]-eigenspace.E0)) # entry in Eigenvector * Boltzmann factor
                    for m1=1:fockstates.norb
                        for m2=1:fockstates.norb
                            nn[1,m1,m2] += weight * fockstates.states[n+1][s][j][2*m1-1] * fockstates.states[n+1][s][j][2*m2-0] 
                            nn[2,m1,m2] += weight * fockstates.states[n+1][s][j][2*m1-1] * fockstates.states[n+1][s][j][2*m2-1] 
                            nn[3,m1,m2] += weight * fockstates.states[n+1][s][j][2*m1-0] * fockstates.states[n+1][s][j][2*m2-0] 
                        end # m2
                    end # m1
                end # j
            end # i
        end # s
    end # n
    return nn/Z
end
