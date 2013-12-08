import scipy, scipy.sparse, scipy.linalg

import cPickle as pickle

import HusimiCactus
import TriangleLattice




def LinearizedEigensystem( generation, lattice = 'cactus', periodic = True ):
    """
    Solve via Chris' simpler linearized deviation matrix
    """
    H = GenerateHamiltonian( generation, lattice, periodic )
    
    eigVals, eigVects = scipy.linalg.eig( H )
    
    return eigVals, eigVects



def GenerateHamiltonian( generation, lattice = 'cactus', periodic = True ):
    """
    Generate the linearized hamiltonian from the adjacency list for the Husimi
    cactus or triangular lattice
    """
    
    if lattice == 'cactus':
        N = HusimiCactus.numNodes( generation )
        adjacencyList = HusimiCactus.AdjacencyList( generation, periodic )
    elif lattice == 'triangle':
        N = TriangleLattice.numNodes( generation )
        adjacencyList = TriangleLattice.AdjacencyList( generation )
    
    #   Trim the adjacency list of redundant pairs
    i = 0
    while (i < len(adjacencyList)):
        x,y = adjacencyList[i]
        if x <= y:
            i += 1
        else:
            tossVal = adjacencyList.pop(i)
    
    H = scipy.zeros( (2*N, 2*N) )
    
    h = scipy.zeros( N )
    J_ij = -1.0
    cosTheta_ij = -.5
    
    for pair in adjacencyList:
        i,j = pair
        
        #   theta_i = h_i * ( - y_i + J_ij * y_j )
        H[ i    , j + N ] += J_ij
        H[ j    , i + N ] += J_ij
        
        #   y_i     = h_i * ( theta_i - J_ij * cos(theta_ij) )
        H[ i + N, j     ] += -1.0 * J_ij * cosTheta_ij
        H[ j + N, i     ] += -1.0 * J_ij * cosTheta_ij
    
        #   h_i = Z_i / 2
        h[i] += .5
        h[j] += .5
    
    #   Distributing the h_i terms
    for i in range( int(N) ):
        H[ i    , i + N ] += -1.0
        H[ i + N, i     ] += 1.0
        
        h_i = h[i]
        H[ i    , : ] *= h_i
        H[ i + N, : ] *= h_i
    
    return H
