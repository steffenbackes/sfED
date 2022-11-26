const CAmatrix    = Matrix #SparseMatrixCSC{FockElement,Int64}   # matrix type for creation/annihilation matrix

"""
    getCmatrix(anni::Int64,subspace1, subspace2)

Generate the matrix for the annihilation operator acting on the subspace with N,S quantum number
we act `subspace1` * c * `subspace2`
`subspace 2` has all basis states with N,S
`subspace 1` has all basis states with N-1,S-dS
now just determine the transformation when acting the annihilation operator `anni` on it
`anni` is the orbital/spin index
"""
function getCmatrix(anni::Int64, 
                    subspace1::Vector{Fockstate}, 
                    subspace2::Vector{Fockstate})::CAmatrix
    dim1 = length(subspace1)
    dim2 = length(subspace2)
    indexI = Int64[]; indexJ = Int64[]; val = Int64[]
    # now annihilate the particle in all basis states in subspace2 and find the corresponding state in subspace1
    for j=1:dim2
        if subspace2[j][anni] == 1   # if there is one particle to annihilate...
            state = copy( subspace2[j] )
            c1sgn = getCsign(anni,state) # count the particles before anni
            state[anni] = 0

            # now find this state in the subspace1
            i = findfirst(map(x-> all(x .== state), subspace1))
            push!(indexI,i); push!(indexJ,j); push!(val,c1sgn)
        end
    end # j
    return collect(sparse(indexI,indexJ,val, dim1,dim2))
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
function getCdagmatrix(crea::Int64, 
                       subspace1::Vector{Fockstate}, 
                       subspace2::Vector{Fockstate})::CAmatrix
    dim1 = length(subspace1)
    dim2 = length(subspace2)
    indexI = Int64[]; indexJ = Int64[]; val = Int64[]

    # now create the particle in all basis states and find the corresponding state in the N,S space
    for j=1:dim2
        if subspace2[j][crea] == 0   # if there is space to create one particle...
            state = copy( subspace2[j] )
            c1sgn = getCsign(crea,state) # count the particles before crea
            state[crea] = 1

            # now find this state in the N,S subspace
            i = findfirst(map(x-> all(x .== state), subspace1))
            push!(indexI,i); push!(indexJ,j); push!(val,c1sgn)
        end
    end # j
    return collect(sparse(indexI,indexJ,val, dim1,dim2))
end
######################################################################


"""
    getCCdaggerMat(flavor,subspace1,subspace2)

Return the creation (flavor>0) or annihilation (flavor<0) matrix for the given flavor and subspaces
"""
function getCCdaggerMat(flavor::Int64, 
                       subspace1::Vector{Fockstate}, 
                       subspace2::Vector{Fockstate})::CAmatrix
    if flavor>0
        return getCdagmatrix(abs(flavor), subspace1, subspace2)
    else
        return getCmatrix(abs(flavor), subspace1, subspace2)
    end
end
