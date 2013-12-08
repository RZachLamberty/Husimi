################################################################################
##                                                                            ##
##  HusimiCactus.py                                                           ##
##                                                                            ##
##      A few functions to calculate facts related to the husimi cactus       ##
##                                                                            ##
################################################################################

import scipy


##########################
##                      ##
##  General Functions   ##
##                      ##
##########################

def numNodes( generation, triangleCentered = True ):
    """
    Given the generation, return the number of nodes for a cactus at that gen
    """
    
    if triangleCentered:
        return 3 * scipy.sum( [2**i for i in range(generation+1) ] )
    else:
        raise ValueError, "Site centered not programmed yet"




######################################
##                                  ##
##  Ternary to index conversions    ##
##                                  ##
######################################

def decToTern( N, length = None ):
    """
    Given an integer, express it in trinary (as a list)
    """
    tern = []
    x = int(N)
    while (x != 0):
        tern.append( x % 3 )
        x /= 3
    
    if (length != None):
        if length < len( tern ):
            raise ValueError, "Suggested length is less than the length of the ternary number"
        else:
            tern += [0]*(length - len(tern))
    
    tern.reverse()
    
    return tuple(tern)


def ternToDec( tern ):
    """
    Convert a ternary number to a decimal number
    """
    return scipy.sum( [ 3**i * tern[-(i+1)] for i in range(len(tern)) ] )


def indexToTernDic( generation ):
    """
    Given a generation, return a dictionary which takes indices from 0 to
    numNodes(generation) - 1 into ternary numbers.
    """
    iToTern = {}
    possTernList = [ decToTern(i, generation+1) for i in range( 3**(generation+1) ) ]
    
    i = 0
    N = numNodes( generation )
    while ( i < N ):
        if orderedTern( possTernList[0] ):
            iToTern[i] = possTernList[0]
            i += 1
            tossTern = possTernList.pop(0)
        else:
            tossTern = possTernList.pop(0)
    
    return iToTern


def ternToIndexDic( it ):
    """
    Just reverse an index-to-ternary dictionary
    """
    ternToI = {}
    
    for i in it:
        ternToI[ it[i] ] = i
    
    return ternToI


def orderedTern( tern ):
    """
    Given a ternary tuple, return true if there are no 2's in the center of the
    ternary number
    """
    for i in range(1, len(tern)-1):
        if (tern[i] == 2) and (tern[i+1] != 2):
            return False
    return True



##########################
##                      ##
##  Adjacency matrix    ##
##                      ##
##########################

def AdjacencyList( generation, periodic = False ):
    """
    Given a generation, return a list of index pairs.  It will be redundant.
    """
    adjList = []
    
    indexToTern = indexToTernDic( generation )
    ternToIndex = ternToIndexDic( indexToTern )
    
    ternList = sorted( ternToIndex.keys() )
    
    for tern in ternList:
        neighborTernList = NeighboringTerns( tern, periodic )
        i = ternToIndex[ tern ]
        for neighborTern in neighborTernList:
            adjList.append( (i, ternToIndex[ neighborTern ]) )
    
    return adjList
    

def NeighboringTerns( tern, periodic = False ):
    """
    given a ternary number at a certain generation, return all the neighbors
    """
    neighbors = []
    listTern = list(tern)
    l = len(tern)
    
    #   Find the first 2 after the start, or the end
    i = 1
    while (i < l) and (tern[i] != 2):
        i += 1
    
    #   Regardless of what i ended up being, we will want to cycle through the
    #   numbers immediately in front of that index slot (if i == l, this is the
    #   final element (on an edge piece) and we'll change it either way; if not
    #   it's the values from the previous generation
    n1 = list(listTern)
    n2 = list(listTern)
    n1[i-1] = (n1[i-1] + 1) % 3
    n2[i-1] = (n2[i-1] + 2) % 3
    neighbors.append( tuple(n1) )
    neighbors.append( tuple(n2) )
    
    #   If i != l, we want to cycle through the changes available at that index
    if i != l:
        n1 = list(listTern)
        n2 = list(listTern)
        n1[i] = (n1[i] + 1) % 3
        n2[i] = (n2[i] + 2) % 3
        neighbors.append( tuple(n1) )
        neighbors.append( tuple(n2) )
    
    #   We decided to include the ability to tie together boundary points. This
    #   is how I've decided to do that
    else:
        if periodic:
            n1 = list(listTern)
            n2 = list(listTern)
            n1[0] = (n1[0] + 1) % 3
            n2[0] = (n2[0] + 2) % 3
            neighbors.append( tuple(n1) )
            neighbors.append( tuple(n2) )
    
    return neighbors
