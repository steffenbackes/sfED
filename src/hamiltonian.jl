const Hamiltonian = SparseMatrixCSC{Complex{Float32},Int64}   # matrix type for Hamiltonian matrix

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
    getTmatrix(pModel)

Return the hoppingmatrix, defined by hand here
TODO: this should be overloaded to accept text file or command line inputs.
"""
function getTmatrix(fockstates::Fockstates,pSimulation::SimulationParameters)
   tmatrix = Array{Float64}(undef,fockstates.norb,fockstates.norb)
   t = pSimulation.t
   if fockstates.norb==2
      # A B dimer with one orbital per site
      tmatrix = -[0 t;
                  t 0]
   elseif fockstates.norb==4
      # A B dimer with two orbitals per site
      tmatrix = -[0 0 t 0;
                  0 0 0 t;
                  t 0 0 0;
                  0 t 0 0]
      # A B
      # C D plaquette with one orbital per site
#      tmatrix = -[0 t t 0;
#                  t 0 0 t;
#                  t 0 0 t;
#                  0 t t 0]
   elseif fockstates.norb==5
      # Use Julian's AIM benchmark system
      # beta=11
      t1 = 0.302313674983807
      t2 = 0.154603780668397
      t3 = 0.302313674983807
      t4 = 0.154603780668397
      #beta = 10
      #t1 = 0.304118914651299
      #t2 = 0.152049589882341
      #t3 = 0.304118914651299
      #t4 = 0.152049589882341
      tmatrix =  [0 t1 t2 t3 t4;
                  t1 0 0 0 0;
                  t2 0 0 0 0;
                  t3 0 0 0 0;
                  t4 0 0 0 0]
   elseif fockstates.norb==6
      # A B C trimer with two orbitals per site
      tmatrix = -[0 0 t 0 0 0;
                  0 0 0 t 0 0;
                  t 0 0 0 t 0;
                  0 t 0 0 0 t;
                  0 0 t 0 0 0;
                  0 0 0 t 0 0]
   elseif fockstates.norb==8
      # A B plaquette with two orbitals per site
      # C D
      tmatrix = -[0 0 t 0 t 0 0 0;
                  0 0 0 t 0 t 0 0;
                  t 0 0 0 0 0 t 0;
                  0 t 0 0 0 0 0 t;
                  t 0 0 0 0 0 t 0;
                  0 t 0 0 0 0 0 t;
                  0 0 t 0 t 0 0 0;
                  0 0 0 t 0 t 0 0]
    else
        println("No t matrix implemented for norb=",fockstates.norb,", set hopping matrix to zero!")
        tmatrix = zeros(Float64,fockstates.norb,fockstates.norb)
   end
   return tmatrix
end

"""
    getUJmatrix(pModel)

Return the Coulomb interaction matrices, defined by hand here
TODO: this should be overloaded to accept text file or command line inputs.
"""
function getUJmatrix(fockstates::Fockstates,pSimulation::SimulationParameters)
   Umatrix = Array{Float64}(undef,fockstates.norb,fockstates.norb)
   Jmatrix = Array{Float64}(undef,fockstates.norb,fockstates.norb)
   U = pSimulation.U
   Up = pSimulation.Up
   J = pSimulation.J

   if (fockstates.norb==1)
      Umatrix[1,1] = U
      Jmatrix[1,1] = 0
   elseif (fockstates.norb==2)
      Umatrix = [U Up;
                 Up U]
      Jmatrix = [0 J;
                 J 0]
   elseif fockstates.norb==4
      # A B dimer
      Umatrix = [U Up 0 0;   # we assume two orbitals per site
                 Up U 0 0;
                 0 0 U Up;
                 0 0 Up U]
      Jmatrix = [0 J 0 0;
                 J 0 0 0;
                 0 0 0 J;
                 0 0 J 0]
   elseif fockstates.norb==5
      # Take julian's AIM benchmark
      Umatrix = [U 0 0 0 0;   # single orbital AIM and 4 bath sites
                 0 0 0 0 0;
                 0 0 0 0 0;
                 0 0 0 0 0;
                 0 0 0 0 0]
      Jmatrix = [0 0 0 0 0;
                 0 0 0 0 0;
                 0 0 0 0 0;
                 0 0 0 0 0;
                 0 0 0 0 0.0]
   elseif fockstates.norb==6
      # A B C trimer
      Umatrix = [U Up 0 0 0 0; # we assume two orbitals per site
                 Up U 0 0 0 0;
                 0 0 U Up 0 0;
                 0 0 Up U 0 0;
                 0 0 0 0 U Up;
                 0 0 0 0 Up U]
      Jmatrix = [0 J 0 0 0 0;
                 J 0 0 0 0 0;
                 0 0 0 J 0 0;
                 0 0 J 0 0 0;
                 0 0 0 0 0 J;
                 0 0 0 0 J 0]
   elseif fockstates.norb==8
      # A B plaquette
      # C D
      Umatrix = [U Up 0 0 0 0 0 0; # we assume two orbitals per site
                 Up U 0 0 0 0 0 0;
                 0 0 U Up 0 0 0 0;
                 0 0 Up U 0 0 0 0;
                 0 0 0 0 U Up 0 0;
                 0 0 0 0 Up U 0 0;
                 0 0 0 0 0 0 U Up;
                 0 0 0 0 0 0 Up U]
      Jmatrix = [0 J 0 0 0 0 0 0;
                 J 0 0 0 0 0 0 0;
                 0 0 0 J 0 0 0 0;
                 0 0 J 0 0 0 0 0;
                 0 0 0 0 0 J 0 0;
                 0 0 0 0 J 0 0 0;
                 0 0 0 0 0 0 0 J;
                 0 0 0 0 0 0 J 0]
    else
        println("No UJ matrix implemented for norb=",fockstates.norb,", set interaction to zero!")
   end
   return Umatrix,Jmatrix
end

"""
    getHopping(istate, jstate, tmatrix)

return the hopping contribution for states `<i| |j>`
"""
function getHopping( istate::Fockstate,
                     jstate::Fockstate,
                     tmatrix::Array{Float64,2},
                     pNumerics::NumericalParameters)
   # we act c^dagger c on the state |j> and check overlap with <i|
   norb = size(tmatrix)[1]
   htmp = 0.0
   for m1=1:norb
      for m2=1:norb
         if ( abs(tmatrix[m1,m2]) > pNumerics.cutoff && m1!=m2 )   # no hopping for same orbitals (this is in eps)
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
                     Umatrix::Array{Float64,2},
                     Jmatrix::Array{Float64,2},
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
                        Jmatrix::Array{Float64,2},
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
                                jstate::Fockstate )
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

      if state==istate
         tmp = a4sgn*a3sgn*c2sgn*c1sgn
      end
   end # if
   return tmp
end


"""
    getHamiltonian(eps,tmatrix,Umatrix,Jmatrix,states)

Constructs hamiltonian from onsite energies `eps`, hoppings `tmatrix`, coulomb interactions `Umatrix` and `Jmatrix`
over precomputed `states`.

"""
function getHamiltonian(eps::Array{Float64,1},tmatrix::Array{Float64,2},
                        Umatrix::Array{Float64,2},Jmatrix::Array{Float64,2}, 
                        mu::Float64,aim::Int64,
                        states::Array{Fockstate,1},
                        pNumerics::NumericalParameters)::Hamiltonian

   # Set up index array for the sparse Hamiltonian Matrix
   HamiltonianElementsI = Int64[]
   HamiltonianElementsJ = Int64[]
   HamiltonianElementsV = ComplexF64[]
   for i=1:length(states)
      Hiitmp = 0.0
      # set the diagonals 
      if aim==1
         Hiitmp += -mu*sum(states[i][1:2])     # chemical potential: this is when doing an AIM one-orbital calculation
      else
         Hiitmp += -mu*sum(states[i])     # chemical potential
      end

      Hiitmp += sum(eps .* states[i] ) # onsite levels
   
      # get Density-Density interaction terms
      Hiitmp += getUdensity(states[i],Umatrix,Jmatrix,pNumerics)
   
      # Now we can set the matrix element
      push!( HamiltonianElementsI, i ); push!( HamiltonianElementsJ, i ); push!( HamiltonianElementsV, Hiitmp )
   
      # now offdiagonals
      for j=i+1:length(states)
         Hijtmp = 0.0
   
         # Hopping terms ###########################################
         Hijtmp += getHopping(states[i], states[j], tmatrix,pNumerics)
   
         # Pair-hopping and spin-flip terms of the Coulomb interaction
         Hijtmp += getUnondensity(states[i], states[j], Jmatrix,pNumerics)
   
         # Now we can set the matrix element
         push!( HamiltonianElementsI, i ); push!( HamiltonianElementsJ, j ); push!( HamiltonianElementsV, Hijtmp )
         push!( HamiltonianElementsI, j ); push!( HamiltonianElementsJ, i ); push!( HamiltonianElementsV, Hijtmp ) # hermitian conjugate
                        
      end # j ket
   end # i bra
   
   # Construct the sparse matrix
   return sparse( HamiltonianElementsI,HamiltonianElementsJ,HamiltonianElementsV )
end

"""
    getEvalEvecs(Hamiltonian)

Diagonalize a given `Hamiltonian`. Eigenvalues will be cast to real, since the Hamiltonian is hermitian.
"""
function getEvalEvecs(hamiltonian::Hamiltonian)::Tuple{Array{Eigenvalue,1}, EigenvectorMatrix }

      HamDense = Hermitian(Matrix(hamiltonian))
      evals = eigvals(HamDense)
      evecs = eigvecs(HamDense)
      return real(evals), evecs
end

"""
    getEvalveclist(eps,tmatrix,Umatrix,Jmatrix,mu,fockstates,pNumerics)

Create the N,S submatrices of the Hamiltonian, solve it and return the Eigenvalues and Vectors in a List
"""
function getEvalveclist(eps::Array{Float64,1},tmatrix::Array{Float64,2},
                        Umatrix::Array{Float64,2},Jmatrix::Array{Float64,2},
                        mu::Float64, aim::Int64,
                        fockstates::Fockstates,
                        pNumerics::NumericalParameters)
   Nmax = fockstates.norb*2

   evallist = [ [] for n=0:Nmax]
   eveclist = [ [] for n=0:Nmax]

   for n=0:Nmax
      for s=1:noSpinConfig(n,Nmax)
         dim = length(fockstates.states[n+1][s])
         print("Constructing Hamiltonian(",dim,"x",dim,"), N=",n,", S=",spinConfig(s,n,Nmax),"... ")

         # now get the Hamiltonian submatrix spanned by all states <i| |j> in the N,S space (sparse matrix)
         hamiltonian = getHamiltonian(eps,tmatrix,Umatrix,Jmatrix,mu,aim,fockstates.states[n+1][s],pNumerics)
         println("Done!")

         # (n==8 && spinConfig(s,n,Nmax)==0) && (writeMatrixGnuplot("hamil.dat",Matrix(hamiltonian));  exit())
   
         print("Diagonalizing Hamiltonian... ")
         evals,evecs = getEvalEvecs(hamiltonian)
         
         push!( evallist[n+1], evals )
         push!( eveclist[n+1], [evecs[:,i] for i=1:length(evals)] )
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
