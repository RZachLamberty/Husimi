import scipy, scipy.sparse, scipy.linalg
import pylab
import time
import os
import cPickle as pickle

import HusimiCactus
import TriangleLattice


##################################
##                              ##
##  Some preliminary functions  ##
##  to make sure the directory  ##
##  is set up correctly         ##
##                              ##
##################################

timingPathTri = 'Triangle/timingDic.pkl'
timingPathCac = 'Cactus/timingDic.pkl'
t = {}

for pathName in [timingPathTri, timingPathCac]:
    if not os.path.isfile( pathName ):
        with open( pathName, 'wb' ) as f:
            pickle.dump( t, f )




######################
##                  ##
##  Get to work!    ##
##                  ##
######################

def HolsteinPrimakoff( generation, lattice = "cactus", periodic = True ):
    """
    Given a generation we want to build the cactus out to construct the matrix
    and return the eigenvalues.  We will do this in the semi-roundabout way
    proposed by Mucciolo, Castro Neto, and Chamon in PRB 69 (214424), so what
    is returned is the true eigenspectrum and an eigenvalue matrix representing
    the rotations of a boguliobov-type transformation.
    """
    if lattice == "cactus":
        H = HusimiHamiltonian( generation, periodic )
    elif lattice == "triangle":
        H = TriangleHamiltonian( generation )
    else:
        raise ValueError, "Options are 'cactus' and 'triangle'"
        
    l = H.shape[0] / 2
    K = H[:l, :l]
    L = H[:l, l:]
    
    squaredDiff = scipy.dot(K, K) - scipy.dot(L, L)
    commutator = scipy.dot(L,K) - scipy.dot(K,L)
    
    if scipy.sum(commutator) == 0.0:
        eigVals, eigVects = scipy.linalg.eigh( squaredDiff )
    else:
        eigVals, eigVects = scipy.linalg.eig( squaredDiff - commutator )
    
    #   The 'real' is not a cheat -- zero values could be negative as a result
    #   of roundoff, this takes that into account.
    eigVals = scipy.real( scipy.sqrt( eigVals ) )
    
    return eigVals, eigVects



######################
##                  ##
##  Husimi Cactus   ##
##                  ##
######################

def HusimiHamiltonian( generation, periodic = True ):
    """
    Construct second-order (in Holstein-Primakoff bosons) contribution to the
    hamiltonian for the heisenberg interaction on the Husimi cactus of
    generation 'generation'
    
    Some conventions for indexes in H:
        
        a_i^d   * b_j^d     -->     [ i    , j + N ]
        a_i^d   * b_j       -->     [ i    , j     ]
        a_i     * b_j^d     -->     [ i + N, j + N ]
        a_i     * b_j       -->     [ i + N, j     ]
    
    """
    N = HusimiCactus.numNodes( generation )
    adjacencyList = HusimiCactus.AdjacencyList( generation, periodic )
    
    #   Trim the adjacency list of redundant pairs
    i = 0
    while (i < len(adjacencyList)):
        x,y = adjacencyList[i]
        if x <= y:
            i += 1
        else:
            tossVal = adjacencyList.pop(i)
    
    H = scipy.zeros( (2*N, 2*N) )
    
    for pair in adjacencyList:
        i,j = pair
        H[ i    , j     ] += 1.0        #   s+ t
        H[ j    , i     ] += 1.0        #   t+ s
        H[ i + N, j + N ] += 1.0        #   s  t+
        H[ j + N, i + N ] += 1.0        #   t  s+
        
        H[ i    , i     ] += 2.0        #   s+ s
        H[ j    , j     ] += 2.0        #   t+ t
        H[ i + N, i + N ] += 2.0        #   s  s+
        H[ j + N, j + N ] += 2.0        #   t  t+
        
        H[ i    , j + N ] += -3.0       #   s+ t+
        H[ j    , i + N ] += -3.0       #   t+ s+
        H[ i + N, j     ] += -3.0       #   s  t
        H[ j + N, i     ] += -3.0       #   s  t

    #H *= S / 8.
    H *= 1 / 8.
    
    return H
    

def HusimiDiagonalizeSparse( H ):
    """
    Given a sparse hamiltonian H, return the eigenvalues and eigenvectors
    """
    return HusimiDiagonalize( H.toarray() )


def HusimiDiagonalize( H ):
    """
    Given a regular hamiltonian H, return the eigenvalues and eigenvectors
    """
    w, v = scipy.linalg.eigh( H )
    return w, v






######################
##                  ##
##  Triangle Latice ##
##                  ##
######################

def TriangleHamiltonian( L, cheat = False ):
    """
    Construct second-order (in Holstein-Primakoff bosons) contribution to the
    hamiltonian for the heisenberg interaction on the triangular lattice of
    'width' L
    
    Some conventions for indexes in H:
        
        a_i^d   * b_j^d     -->     [ i    , j + N ]
        a_i^d   * b_j       -->     [ i    , j     ]
        a_i     * b_j^d     -->     [ i + N, j + N ]
        a_i     * b_j       -->     [ i + N, j     ]
    
    """
    N = TriangleLattice.numNodes( L )
    adjacencyList = TriangleLattice.AdjacencyList( L )
    
    #   Trim the adjacency list of redundant pairs
    i = 0
    while (i < len(adjacencyList)):
        x,y = adjacencyList[i]
        if x <= y:
            i += 1
        else:
            tossVal = adjacencyList.pop(i)
    
    H = scipy.zeros( (2*N, 2*N) )
    
    for pair in adjacencyList:
        i,j = pair
        H[ i    , j     ] += 1.0        #   s+ t
        H[ j    , i     ] += 1.0        #   t+ s
        H[ i + N, j + N ] += 1.0        #   s  t+
        H[ j + N, i + N ] += 1.0        #   t  s+
        
        H[ i    , i     ] += 2.0        #   s+ s
        H[ j    , j     ] += 2.0        #   t+ t
        H[ i + N, i + N ] += 2.0        #   s  s+
        H[ j + N, j + N ] += 2.0        #   t  t+
        
        H[ i    , j + N ] += -3.0       #   s+ t+
        H[ j    , i + N ] += -3.0       #   t+ s+
        H[ i + N, j     ] += -3.0       #   s  t
        H[ j + N, i     ] += -3.0       #   s  t

    #H *= S / 8.
    H *= 1 / 8.
    
    if cheat:
        mV = max( H.diagonal() )
        for i in range(len(H)):
            H[i,i] = mV
    
    return H
    

def TriangleDiagonalizeSparse( H ):
    """
    Given a sparse hamiltonian H, return the eigenvalues and eigenvectors
    """
    return TriangleDiagonalize( H.toarray() )


def TriangleDiagonalize( H ):
    """
    Given a regular hamiltonian H, return the eigenvalues and eigenvectors
    """
    w, v = scipy.linalg.eigh( H )
    return w, v



######################
##                  ##
##  Timing stuff    ##
##                  ##
######################

def timeTrial( trialDic, generation ):
    """
    Record how long it takes to calculate the eigenvectors and values for a
    given generation, and save it to a trial dictionary
    """
    if generation in trialDic:
        pass
    else:
        startTime = time.time()
        eigVals, eigVects = HolsteinPrimakoff( generation )
        endTime = time.time()
        
        trialDic[ generation ] = endTime - startTime




##################################
##                              ##
##  Doing the calculations      ##
##                              ##
##################################

def RecordThroughGen( generation, lattice = "cactus", startGen = 1, periodic = True ):
    """
    Record the eigenvalues and eigenvectors up through a generation
    """
    for gen in range( startGen, generation + 1 ):
        print "\tGeneration " + str( gen ).zfill(2) + "..."
        
        if lattice == "cactus":
            prefix = 'Cactus/'
        elif lattice == "triangle":
            prefix = 'Triangle/'
        eigValFile = prefix + 'eigVal_gen' + str(gen).zfill(2) + '.pkl'
        eigVecFile = prefix + 'eigVec_gen' + str(gen).zfill(2) + '.pkl'
        
        if periodic:
            eigValFile = eigValFile.replace( '.pkl', '_periodic.pkl' )
            eigVecFile = eigVecFile.replace( '.pkl', '_periodic.pkl' )
        
        if os.path.isfile( eigValFile ) and os.path.isfile( eigVecFile ):
            pass
        
        else:
            startTime = time.time()
            eigVals, eigVecs = HolsteinPrimakoff( gen, lattice, periodic )
            endTime = time.time()
            
            #   pickle the eigenvalues and eigenvectors
            with open( eigValFile, 'wb' ) as f:
                pickle.dump( eigVals, f )
            with open( eigVecFile, 'wb' ) as f:
                pickle.dump( eigVecs, f )
            
            #   Record timing info
            timingFile = prefix + 'timingDic.pkl'
            with open( timingFile, 'rb' ) as f:
                t = pickle.load( f )
                t[gen] = endTime - startTime
            with open( timingFile, 'wb' ) as f:
                pickle.dump(t, f)


def LoadEigenValues(gen = 'max', lattice = 'cactus', periodic = True):
    
    if lattice == 'cactus':
        prefix = "Cactus/"
    elif lattice == 'triangle':
        prefix = "Triangle/"
    
    if gen == "max":
        gen = MaxGenRecorded( lattice )
        
    evPath = prefix + 'eigVal_gen' + str(gen).zfill(2)
    
    if periodic:
        evPath += '_periodic.pkl'
    else:
        evPath += '.pkl'
        
    if os.path.isfile(evPath):
        with open( evPath, 'rb' ) as f:
            evs = pickle.load( f )
        return evs
    else:
        raise ValueError, "Generation hasn't been calculated yet"


def LoadEigenVectors(gen = 'max', lattice = 'cactus', periodic = True):
    
    if lattice == 'cactus':
        prefix = "Cactus/"
    elif lattice == 'triangle':
        prefix = "Triangle/"
    
    if gen == "max":
        gen = MaxGenRecorded( lattice )
        
    evPath = prefix + 'eigVec_gen' + str(gen).zfill(2)
    
    if periodic:
        evPath += '_periodic.pkl'
    else:
        evPath += '.pkl'
        
    if os.path.isfile(evPath):
        with open( evPath, 'rb' ) as f:
            evs = pickle.load( f )
        return evs
    else:
        raise ValueError, "Generation hasn't been calculated yet"


def MaxGenRecorded( lattice = 'cactus', periodic = True ):
    
    if lattice == 'cactus':
        prefix = "Cactus/"
    elif lattice == 'triangle':
        prefix = "Triangle/"
    
    fList = os.listdir( prefix )
    fList = [int(x[10:12]) for x in fList if x.startswith('eigVal') and x.endswith('periodic.pkl')]
    
    return max(fList)


def LoadTimingInfo(lattice = 'cactus'):
    
    if lattice == 'cactus':
        prefix = "Cactus/"
    elif lattice == 'triangle':
        prefix = "Triangle/"
        
    tPath = prefix + 'timingDic.pkl'
    if os.path.isfile(tPath):
        with open( tPath, 'rb' ) as f:
            t = pickle.load( f )
        return t
    else:
        raise ValueError, "No timing info exists yet"


def PlotTimingScaling( lattice = 'cactus'):
    
    if lattice == 'cactus':
        prefix = "Cactus/"
    elif lattice == 'triangle':
        prefix = "Triangle/"
        
    t = LoadTimingInfo( lattice )
    x = scipy.array( sorted( t.keys() ) )
    y = scipy.array( [t[i] for i in x] )
    
    f1 = pylab.figure()
    p1 = f1.add_subplot(111)
    p1.plot( x, y / 60. )
    
    f2 = pylab.figure()
    p2 = f2.add_subplot(111)
    p2.plot( x, y / 60. )
    p2.set_yscale( 'log' )
