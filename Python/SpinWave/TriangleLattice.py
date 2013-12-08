################################################################################
##                                                                            ##
##  TriangleLattice.py                                                        ##
##                                                                            ##
##      A few functions to calculate facts related to the triangle lattice    ##
##                                                                            ##
################################################################################

import scipy


##########################
##                      ##
##  General Functions   ##
##                      ##
##########################

def numNodes( L ):
    """
    Given the side length L, return the number of nodes for a triangular lattice
    of 'width' L (i.e. a square L x L grid)
    """
    return L**2.


##########################
##                      ##
##  Adjacency matrix    ##
##                      ##
##########################

def AdjacencyList( L ):
    """
    Given a 'width' L, return a list of index pairs.  It will be redundant.
    """
    adjList = []
    
    for i in range(L**2):
        neighborList = NeighboringNodes( i, L )
        for neighbor in neighborList:
            adjList.append( (i, neighbor) )
    
    return adjList
    

def NeighboringNodes( i, L ):
    """
    given an index i, return all the neighbors in an L x L triangular lattice
    """
    neighbors = []
    
    ix = i % L
    iy = i / int(L)
    
    #   Up
    neighbors.append( ix                +   L * ( ( iy + 1 ) % L )  )
    
    #   Down
    neighbors.append( ix                +   L * ( ( iy - 1 ) % L )  )
    
    #   Right
    neighbors.append( ( ix + 1 ) % L    +   L * ( iy )              )
    
    #   Left
    neighbors.append( ( ix - 1 ) % L    +   L * ( iy )              )
    
    #   Up-right
    neighbors.append( ( ix + 1 ) % L    +   L * ( ( iy + 1 ) % L )  )
    
    #   Down-left
    neighbors.append( ( ix - 1 ) % L    +   L * ( ( iy - 1 ) % L )  )
    
    return neighbors
