################################################################################
##                                                                            ##
##  makeHam.py                                                                ##
##                                                                            ##
##      This file will generate a hamiltonian for the Heisenberg interaction  ##
##      between spin 1/2 elements on a Husimi cactus.  At first it will be an ##
##      exact diagonalization routine for generation 1 (4 triangles)          ##
##                                                                            ##
################################################################################

import scipy, scipy.sparse, scipy.linalg


##########################
##                      ##
##  Site basis only     ##
##                      ##
##########################


def exactEig( gen, Jf, Ja):
    """
    Given the generation and values of Jf and Ja, return the eigenvalues and
    eigenvectors of the Heisenberg interaction on that generation's Husimi
    cactus
    """
    H = makeH( gen, Jf, Ja )
    return scipy.linalg.eigh( H.toarray() )


def makeH( gen, Jf, Ja ):
    """
    return a sparse array (dictionary of (conf_i, conf_j) : h_ij pairs) that
    represents the Heisenberg hamiltonian on the husimi cactus
    """

    numNodes = numberOfNodes( gen )

    H = scipy.sparse.dok_matrix( (2**numNodes, 2**numNodes), dtype = scipy.float32 )
    
    for i in range( 2**numNodes ):
        
        repString = binaryRep( i, gen )
        
        triAdj = triangleAdjacency( gen )
        edgeAdj = edgeAdjacency( gen )
        
        triHam( H, Ja, i, repString, triAdj )
        
        edgeHam( H, Jf, i, repString, edgeAdj )
    
    return H


def numberOfNodes( gen ):
    """
    Given the generation, return the number of nodes in the Husimi cactus at
    that generation (gen 0 is one triangle, gen 1 is 4, etc)
    """
    return int( scipy.sum( [ 3.**i for i in range( 1, gen + 2 ) ] ) )


def binaryRep( i, gen ):
    """
    given an integer and a generation, return the representation of that
    configuration i in the site basis, as a string of 0's and 1's
    
    Basically, at gen we should have numberOfNodes points in our binary rep of
    i.
    """
    
    length = numberOfNodes( gen )
    b = scipy.binary_repr( i, length )
    
    return b
    

def triangleAdjacency( gen ):
    """
    return a list of triangle adjacency index tuples for generation gen
    """
    
    numTri = int( numberOfNodes( gen ) / 3. )
    
    return [ ( 3*i+j, 3*i+((j+1)%3) ) for j in range(3) for i in range(numTri) ]


def edgeAdjacency( gen ):
    """
    return a list of edge adjacency index tuples for generation gen
    """
    if gen == 0:
        return []
    elif gen == 1:
        return [(0,5), (1,8), (2,11)]
    else:
        raise ValueError, "Hasn't been programmed yet!"


def triHam( H, Ja, repNum, repString, adjacencyList ):
    """
    Given a hamiltonian to update, a representative number and string, and a
    list of neighbors on triangles, do your stuff!
    """
    
    flipVal = .5 * Ja
    opVal = -.25 * Ja
    sameVal = .25 * Ja
    
    for (i,j) in adjacencyList:
        si,sj = repString[i], repString[j]
        
        if si != sj:
            #   Raising / Lowering terms will flip the two and spit out 1/2 Ja
            #   We have to do some wonky stuff because you can't do element 
            #   assignment with strings
            flipString = [c for c in repString]
            flipString[i] = sj
            flipString[j] = si
            flipString = ''.join( flipString )
            flipNum = int( flipString, base = 2 )
            
            H[(repNum, flipNum)] += flipVal
            
            #   Opposite spins will get a -1/4 from the s_z terms
            H[(repNum, repNum)] += opVal
            
        else:
            #   Same spins will get a +1/4 from the s_z terms
            H[(repNum, repNum)] += sameVal
    

def edgeHam( H, Jf, repNum, repString, adjacencyList ):
    """
    Given a hamiltonian to update, a representative number and string, and a
    list of neighbors on connecting edges, do your stuff!
    """
    
    flipVal = -.5 * Jf
    opVal = .25 * Jf
    sameVal = -.25 * Jf
    
    #   Cycle throught the adjacencies
    for (i,j) in adjacencyList:
        
        si,sj = repString[i], repString[j]
        
        if si != sj:
            #   Raising / Lowering terms will flip the two and spit out -1/2 Jf
            #   We have to do some wonky stuff because you can't do element 
            #   assignment with strings
            flipString = [c for c in repString]
            flipString[i] = sj
            flipString[j] = si
            flipString = ''.join( flipString )
            flipNum = int( flipString, base = 2 )
            
            H[(repNum, flipNum)] += flipVal
            
            #   Opposite spins will get a 1/4 from the s_z terms
            H[(repNum, repNum)] += opVal
            
        else:
            #   Same spins will get a -1/4 from the s_z terms
            H[(repNum, repNum)] += sameVal



##############################
##                          ##
##  Local Basis instead     ##
##                          ##
##############################


def simpleH(Ja):
    """
    The Hamiltonian of just one triangle
    """
    
    Ja = float(Ja)
    return scipy.array( [                                                   \
    [ 3*Ja/4,   Ja/2,   Ja/2,   Ja/2,   0,      0,      0,      0       ],  \
    [ Ja/2,     -Ja/4,  0,      0,      0,      Ja/2,   Ja/2,   0       ],  \
    [ Ja/2,     0,      -Ja/4,  0,      Ja/2,   0,      Ja/2,   0       ],  \
    [ Ja/2,     0,      0,      -Ja/4,  Ja/2,   Ja/2,   0,      0       ],  \
    [ 0,        0,      Ja/2,   Ja/2,   -Ja/4,  0,      0,      Ja/2    ],  \
    [ 0,        Ja/2,   0,      Ja/2,   0,      -Ja/4,  0,      Ja/2    ],  \
    [ 0,        Ja/2,   Ja/2,   0,      0,      0,      -Ja/4,  Ja/2    ],  \
    [ 0,        0,      0,      0,      Ja/2,   Ja/2,   Ja/2,   3*Ja/4  ]   \
    ])


def truncBasisH(Ja, num):
    """
    Given a value for Ja, return the lowest num energy states and eigenvectors.
    w[i] <----> v[:,i]
    """
    w, v = scipy.linalg.eigh( simpleH( Ja ) )
    minNum = min( num, len( w ) )
    return w[:minNum], v[:,:minNum]


def sPlusAndMinusAndZ( v ):
    """
    given a matrix of eigenvectors v[:,i] in the site basis, calculate s_i^+ and
    s_i^-
    """
    numN = v.shape[1]
    sp = scipy.zeros((3, numN, numN))
    sm = scipy.zeros((3, numN, numN))
    sz = scipy.zeros((3, numN, numN))
    
    spV = scipy.zeros((3,) + v.shape)
    smV = scipy.zeros((3,) + v.shape)
    szV = scipy.zeros((3,) + v.shape)
    
    #   Calculate Si+ | v_i > and Si- | v_i > and Siz |v_j>
    for j in range(8):
        b = '{0:03b}'.format(j)
        
        #   Splus and Sminus and Sz
        
        #   Site 0
        if b[0] == '0':
            bNew = '1' + b[1:]
            jNew = int( bNew, 2 )
            spV[0, jNew, :] = v[jNew, :]
            szV[0, j, :] = -.5 * v[j, :]
        else:
            bNew = '0' + b[1:]
            jNew = int( bNew, 2 )
            smV[0, jNew, :] = v[jNew, :]
            szV[0, j, :] = .5 * v[j, :]
        
        #   Site 1
        if b[1] == '0':
            bNew = b[0] + '1' + b[2]
            jNew = int( bNew, 2 )
            spV[1, jNew, :] = v[jNew, :]
            szV[1, j, :] = -.5 * v[j, :]
        else:
            bNew = b[0] + '1' + b[2]
            jNew = int( bNew, 2 )
            smV[1, jNew, :] = v[jNew, :]
            szV[1, j, :] = .5 * v[j, :]
            
        #   Site 2
        if b[2] == '0':
            bNew = b[:2] + '1'
            jNew = int( bNew, 2 )
            spV[2, jNew, :] = v[jNew, :]
            szV[2, j, :] = -.5 * v[j, :]
        else:
            bNew = b[:2] + '0'
            jNew = int( bNew, 2 )
            smV[2, jNew, :] = v[jNew, :]
            szV[2, j, :] = .5 * v[j, :]
    
    #   calculate the products
    for i in range(3):
        for n in range(numN):
            for np in range(n, numN):
                xp = scipy.dot( v[:,n], spV[i, :, np] )
                sp[i, n, np] = xp
                sp[i, np, n] = xp
                xm = scipy.dot( v[:,n], smV[i, :, np] )
                sm[i, n, np] = xm
                sm[i, np, n] = xm
                xz = scipy.dot( v[:,n], szV[i, :, np] )
                sz[i, n, np] = xz
                sz[i, np, n] = xz
    
    return sp, sm, sz


def glueEmH( Ja, Jf, truncNum = scipy.inf ):
    """
    Given Ja, Jf, and the number of states we wish to keep at each step, glue
    together four copies with bonds betwixt, and return the Hamiltonian of the
    combination
    """
    w, v = truncBasisH( Ja, truncNum )
    sPlus, sMinus, sZ = sPlusAndMinusAndZ( v )
    
    H1 = scipy.zeros( ( len(w)**4, len(w)**4 ) )
    
    for n in range( len(w)**4 ):
        #   Diagonal previous generation contributions
        o = oct(n)[-4:].zfill(4)
        o = [int(char) for char in o]
        o_A, o_B, o_C, o_D = o
        
        H1[n, n] += scipy.sum( [ w[ i ] for i in o ] )
        
        #   Edge terms
        for np in range( n, len(w)**4 ):
            op = oct(np)[-4:].zfill(4)
            op = [int(char) for char in op]
            op_A, op_B, op_C, op_D = op
            
            x = 0.
            if ( (o_B == op_B) and (o_C == op_C) ):
                x += -Jf * ( .5 * ( sPlus[0][o_A, op_A] * sMinus[0][o_D, op_D] + sMinus[0][o_A, op_A] * sPlus[0][o_D,op_D] ) + sZ[0][o_A, op_A] * sZ[0][o_D, op_D] )
            if ( (o_C == op_C) and (o_A == op_A) ):
                x += -Jf * ( .5 * ( sPlus[1][o_B, op_B] * sMinus[1][o_D, op_D] + sMinus[1][o_B, op_B] * sPlus[1][o_D,op_D] ) + sZ[1][o_B, op_B] * sZ[1][o_D, op_D] )
            if ( (o_A == op_A) and (o_B == op_B) ):
                x += -Jf * ( .5 * ( sPlus[2][o_C, op_C] * sMinus[2][o_D, op_D] + sMinus[2][o_C, op_C] * sPlus[1][o_D,op_D] ) + sZ[1][o_C, op_C] * sZ[2][o_D, op_D] )
            
            H1[n, np] = x
            H1[np, n] = x
    
    return H1


def eigenCheat( Ja, Jf, truncNum = scipy.inf ):
    """
    Given the above, calculate the eigenvalues for the "first generation" cactus
    using the DMRG-esque gluing method
    """
    H = glueEmH( Ja, Jf, truncNum )
    
    return scipy.linalg.eigh( H )
