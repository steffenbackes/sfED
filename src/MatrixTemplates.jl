"""
    getTmatrix(pModel)

Return the hoppingmatrix, defined by hand here
TODO: this should be overloaded to accept text file or command line inputs.
"""
function getTmatrix(fockstates::Fockstates,pSimulation::ModelParameters)
    t = pSimulation.t
    tmatrix = if fockstates.norb==2
        # A B dimer with one orbital per site
        @SMatrix [0 -t;
                  -t 0]
    elseif fockstates.norb==4
        # A B dimer with two orbitals per site
        @SMatrix [ 0  0 -t  0;
                   0  0  0 -t;
                  -t  0  0  0;
                   0 -t  0  0]
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
        @SMatrix [0 t1 t2 t3 t4;
                  t1 0  0  0  0;
                  t2 0  0  0  0;
                  t3 0  0  0  0;
                  t4 0  0  0  0]
    elseif fockstates.norb==6
        # A B C trimer with two orbitals per site
        @SMatrix  [ 0  0 -t  0  0  0;
                    0  0  0 -t  0  0;
                   -t  0  0  0 -t  0;
                    0 -t  0  0  0 -t;
                    0  0 -t  0  0  0;
                    0  0  0 -t  0  0]
    elseif fockstates.norb==8
        # A B plaquette with two orbitals per site
        # C D
        @SMatrix [ 0  0 -t  0 -t  0  0  0;
                   0  0  0 -t  0 -t  0  0;
                  -t  0  0  0  0  0 -t  0;
                   0 -t  0  0  0  0  0 -t;
                  -t  0  0  0  0  0 -t  0;
                   0 -t  0  0  0  0  0 -t;
                   0  0 -t  0 -t  0  0  0;
                   0  0  0 -t  0 -t  0  0]
    else
        println("No t matrix implemented for norb=",fockstates.norb,", set hopping matrix to zero!")
        @SMatrix zeros(Float64,fockstates.norb,fockstates.norb)
    end
    return tmatrix
end

"""
    getUJmatrix(pModel)

Return the Coulomb interaction matrices, defined by hand here
TODO: this should be overloaded to accept text file or command line inputs.
"""
function getUJmatrix(fockstates::Fockstates,pSimulation::ModelParameters)
    U::Float64   = pSimulation.U
    Up::Float64  = pSimulation.Up
    J::Float64   = pSimulation.J

    Umatrix, Jmatrix = if (fockstates.norb==1)
        SMatrix{1,1}([U;;]), SMatrix{1,1}([0.0;;])
    elseif (fockstates.norb==2)
        SMatrix{2,2}([U Up;
                      Up U]),
        SMatrix{2,2}([0 J;
                      J 0])
    elseif fockstates.norb==4
        # A B dimer
        SMatrix{4,4}([U Up 0 0;   # we assume two orbitals per site
                      Up U 0 0;
                      0 0 U Up;
                      0 0 Up U]),
        SMatrix{4,4}([0 J 0 0;
                      J 0 0 0;
                      0 0 0 J;
                      0 0 J 0])
    elseif fockstates.norb==5
        # Take julian's AIM benchmark
        SMatrix{5,5}([U 0 0 0 0;   # single orbital AIM and 4 bath sites
                      0 0 0 0 0;
                      0 0 0 0 0;
                      0 0 0 0 0;
                      0 0 0 0 0]),
        SMatrix{5,5}([0 0 0 0 0;
                      0 0 0 0 0;
                      0 0 0 0 0;
                      0 0 0 0 0;
                      0 0 0 0 0.0])
    elseif fockstates.norb==6
        # A B C trimer
        SMatrix{6,6}([U Up 0 0 0 0; # we assume two orbitals per site
                      Up U 0 0 0 0;
                      0 0 U Up 0 0;
                      0 0 Up U 0 0;
                      0 0 0 0 U Up;
                      0 0 0 0 Up U]),
        SMatrix{6,6}([0 J 0 0 0 0;
                      J 0 0 0 0 0;
                      0 0 0 J 0 0;
                      0 0 J 0 0 0;
                      0 0 0 0 0 J;
                      0 0 0 0 J 0])
    elseif fockstates.norb==8
        # A B plaquette
        # C D
        SMatrix{8,8}([U Up 0 0 0 0 0 0; # we assume two orbitals per site
                      Up U 0 0 0 0 0 0;
                      0 0 U Up 0 0 0 0;
                      0 0 Up U 0 0 0 0;
                      0 0 0 0 U Up 0 0;
                      0 0 0 0 Up U 0 0;
                      0 0 0 0 0 0 U Up;
                      0 0 0 0 0 0 Up U]),
        SMatrix{8,8}([0 J 0 0 0 0 0 0;
                      J 0 0 0 0 0 0 0;
                      0 0 0 J 0 0 0 0;
                      0 0 J 0 0 0 0 0;
                      0 0 0 0 0 J 0 0;
                      0 0 0 0 J 0 0 0;
                      0 0 0 0 0 0 0 J;
                      0 0 0 0 0 0 J 0])
    else
        no = fockstates.norb
        println("No UJ matrix implemented for norb=",no,", set interaction to zero!")
        SMatrix{no,no}(zeros(no,no)),SMatrix{no,no}(zeros(no,no))
    end
    return Umatrix,Jmatrix
end
