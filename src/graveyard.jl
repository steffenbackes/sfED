function getHamiltonian(eps::Vector{Float64},tmatrix::SMatrix,
                        Umatrix::SMatrix,Jmatrix::SMatrix, 
                        mu::Float64,aim::Int64,
                        states::Vector{Fockstate},
                        pNumerics::NumericalParameters)::Hamiltonian

    Hsize = length(states)
    #println("states: ", states)
    H = if Hsize < 2
        #Hi = MMatrix{Hsize,Hsize}(zeros(ComplexF64, Hsize, Hsize))
        Hi = Matrix{ComplexF64}(undef, Hsize, Hsize)
        for i=1:Hsize
            # set the diagonals (chemical potential differs between AIM and Hubbard model)
            Hi[i,i]  = (aim == 1) ? -mu*sum(states[i][1:2]) : -mu*sum(states[i])
            Hi[i,i] += sum(eps .* states[i])           # onsite levels

            # get Density-Density interaction terms
            Hi[i,i] += getUdensity(states[i], Umatrix, Jmatrix, pNumerics)

            # now offdiagonals
            for j=i+1:Hsize
                # Hopping terms ###########################################
                Hijtmp  = getHopping(states[i], states[j], tmatrix,pNumerics)
                # Pair-hopping and spin-flip terms of the Coulomb interaction
                Hijtmp += getUnondensity(states[i], states[j], Jmatrix,pNumerics)

                Hi[i,j] = Hijtmp 
                Hi[j,i] = conj(Hijtmp)
            end # j ket
        end # i bra
        Hi
    else # large martrix, use sparse!!
        # Set up index array for the sparse Hamiltonian Matrix
        HamiltonianElementsI = Stack{Int64}()
        HamiltonianElementsJ = Stack{Int64}()
        HamiltonianElementsV = Stack{ComplexF64}()
        for i=1:Hsize
            # set the diagonals  (chemical potential)
            Hiitmp = (aim==1) ? -mu*sum(states[i][1:2]) : -mu*sum(states[i])
            Hiitmp += sum(eps .* states[i])           # onsite levels
            # get Density-Density interaction terms
            Hiitmp += getUdensity(states[i],Umatrix,Jmatrix,pNumerics)
            # Now we can set the matrix element
            if Hiitmp != 0.0
                push!(HamiltonianElementsI, i); push!(HamiltonianElementsJ, i); push!(HamiltonianElementsV, Hiitmp)
            end

            # now offdiagonals
            for j=i+1:Hsize
                # Hopping terms ###########################################
                Hijtmp  = getHopping(states[i], states[j], tmatrix,pNumerics)
                # Pair-hopping and spin-flip terms of the Coulomb interaction
                Hijtmp += getUnondensity(states[i], states[j], Jmatrix,pNumerics)

                # Now we can set the matrix element
                if Hiitmp != 0.0
                    push!(HamiltonianElementsI, i); push!(HamiltonianElementsJ, j); push!(HamiltonianElementsV, Hiitmp)
                    push!(HamiltonianElementsI, j); push!(HamiltonianElementsJ, i); push!(HamiltonianElementsV, conj(Hiitmp))
                end
            end # j ket
        end # i bra
        H = sparse(collect(HamiltonianElementsI), collect(HamiltonianElementsJ), collect(HamiltonianElementsV))
        println(nnz(H), " vs ", prod(size(H)))
        H
    end

    # Construct the sparse matrix
    println("=======================")
    println(tmatrix)
    println(Umatrix)
    println(Jmatrix)
    println(Matrix(H))
    println("=======================")
    return H
end
