################################################################################
##                                                                            ##
##  BowTieCheck.py.                                                           ##
##                                                                            ##
##      We are going to have preliminary results on the husimi cactus soon,   ##
##      and we hope to have some values to check against.                     ##
##                                                                            ##
################################################################################

import scipy, scipy.linalg
import pylab
import os
import cPickle as pickle


def BowTieEigenSys( S = .5, Jx = 1.0, Jz = 1.0, d = 1.0 ):
    """
    Generate the hamiltonian
    
        H = \sum_{<i,j>} ( Jx * (1/2)[ S+_i S-_j + S-_i S+_j ] + Jz Sz_i Sz_j )
            
            + \sum_{i} ( d * Sz_i ^ 2 )
    
    on a single bowtie, and then diagonalize it.
    """
    H = BowTieH( S, Jx, Jz, d )
    
    if scipy.sum(scipy.conj(H.transpose()) - H) != 0:
        raise ValueError, "Doesn't look to be Hermitian -- check it!"
    else:
        return scipy.linalg.eigh( H )


def BowTieH( S = .5, Jx = 1.0, Jz = 1.0, d = 1.0 ):
    """
    The Hamiltonian for a single bowtie, (2*S + 1)^5 states
    """
    
    A = BowTieAdjacencyDic()
    
    N = (2*S + 1)**5
    
    sPlusMinus, sMinusPlus, sZZ, sZ2 = HParts( S, A )
    
    H = Jx * .5 * ( sPlusMinus + sMinusPlus ) + Jz * sZZ + d * sZ2
    
    return H


def BowTieEigenSys_paramSet( paramSet ):
    """
    Generate the hamiltonian
    
        H = \sum_{<i,j>} ( Jx * (1/2)[ S+_i S-_j + S-_i S+_j ] + Jz Sz_i Sz_j )
            
            + \sum_{i} ( d * Sz_i ^ 2 )
    
    on a single bowtie, and then diagonalize it, for all of the parameters
    provided in paramSet
    """
    for (H, params) in BowTieH_paramSet( paramSet ):
        if scipy.sum( scipy.conj( H.transpose() ) - H) != 0:
            print "Doesn't look to be Hermitian -- might want to check it!"
            w, v = scipy.linalg.eig( H )
            yield (w, v, params)
        else:
            w, v = scipy.linalg.eigh( H )
            yield (w, v, params)


def BowTieH_paramSet( paramSet ):
    """
    The Hamiltonian for a single bowtie, (2*S + 1)^5 states, for a set of d
    values -- this will be a generator equation
    """
    
    A = BowTieAdjacencyDic()
    
    #   Break the paramater set into sub-sets based on the spin value.
    paramSetDic = {}
    paramSet = scipy.array( paramSet )
    
    for row in paramSet:
        sNow = row[0]
        if sNow not in paramSetDic:
            paramSetDic[ sNow ] = []
        
        paramSetDic[ sNow ].append( row )
    
    for S in paramSetDic:
        paramSetDic[S] = scipy.array( paramSetDic[S] )
    
        N = (2*S + 1)**5
        
        sPlusMinus, sMinusPlus, sZZ, sZ2 = HParts( S, A )
        
        for params in paramSetDic[S]:
            S, Jx, Jz, d = params
            H = Jx * .5 * ( sPlusMinus + sMinusPlus ) + Jz * sZZ + d * sZ2
            yield (H, params)


def HParts( S, A ):
    """
    Return the individual parts of the Hamiltonian
    """
    
    N = int(round((2*S + 1)))**5
    sPlusMinus = scipy.zeros( (N, N) )
    sMinusPlus = scipy.zeros( (N, N) )
    sZZ = scipy.zeros( (N, N) )
    sZ2 = scipy.zeros( (N, N) )
    
    #   Configurations are in (2S+1)-nary
    for index in range( N ):
        
        config = IndexToVector( S, index )
        
        #   Cycle through by bond for the first three terms of the hamiltonian
        for i in A:
            for j in A[i]:
                
                #   S+_i S-_j
                connecIndex, fac = SPlusMinus( S, config, i, j )
                if connecIndex != None:
                    sPlusMinus[connecIndex, index] += fac
                
                #   S-_i S+_j
                connecIndex, fac = SMinusPlus( S, config, i, j )
                if connecIndex != None:
                    sMinusPlus[connecIndex, index] += fac
                
                #   Sz_i Sz_j
                sZZ[index, index] += SZZ( S, config, i, j )
            
        #   Sz_i ^ 2
        for i in range( 5 ):
            sZ2[index, index] += SZ2( S, config, i )
    
    return sPlusMinus, sMinusPlus, sZZ, sZ2


def SPlusMinus( S, config, i, j ):
    """
    Return the index of the configuration connected to that given by S+_i S-_j.
    If none exists, return None
    """
    sI = config[i]
    sJ = config[j]
    
    if (sI == int( round(2*S) ) ) or (sJ == 0):
        return None, None
    else:
        copyConfig = list(config)
        copyConfig[i] += 1
        copyConfig[j] -= 1
        
    #   Calculating the Sz factors ( S+- |S Sz > = sqrt( (S -+ Sz)(S +- Sz + 1) ) | S Sz+-1 >
    #fac_i = scipy.sqrt( ( S - (sI - S) ) * ( S + (sI - S) + 1.) )
    #fac_j = scipy.sqrt( ( S + (sI - S) ) * ( S - (sI - S) + 1.) )
    fac_i = scipy.sqrt( ( 2*S - sI ) * ( sI       + 1.) )
    fac_j = scipy.sqrt( ( sJ       ) * ( 2*S - sJ + 1.) )
    
    return (VectorToIndex( S, copyConfig ), fac_i * fac_j)


def SMinusPlus( S, config, i, j ):
    """
    Return the index of the configuration connected to that given by S-_i S+_j.
    If none exists, return None
    """
    sI = config[i]
    sJ = config[j]
    
    if (sI == 0) or (sJ == int( round(2*S) )):
        return None, None
    else:
        copyConfig = list(config)
        copyConfig[i] -= 1
        copyConfig[j] += 1
    
    #   Calculating the Sz factors ( S+- |S Sz > = sqrt( (S -+ Sz)(S +- Sz + 1) ) | S Sz+-1 >
    #fac_i = scipy.sqrt( ( S + (sI - S) ) * ( S - (sI - S) + 1.) )
    #fac_j = scipy.sqrt( ( S - (sI - S) ) * ( S + (sI - S) + 1.) )
    fac_i = scipy.sqrt( ( sI       ) * ( 2*S - sI + 1.) )
    fac_j = scipy.sqrt( ( 2*S - sJ ) * ( sJ       + 1.) )
    
    return (VectorToIndex( S, copyConfig ), fac_i * fac_j)
    

def SZZ( S, config, i, j ):
    """
    Return the value of Sz_i * Sz_j
    """
    return (config[i] - S) * (config[j] - S)


def SZ2( S, config, i ):
    """
    Return Sz_i ^ 2
    """
    return (config[i] - S)**2.


def BowTieAdjacencyDic():
    """
    Return the adjacency dictionary of the bow tie with indices given by
    
        1 ----- 2
         \     /
          \   /
           \ /
            0
           / \
          /   \
         /     \
        3 ----- 4
    
    To make it non-redundant we will require that
        
        i ---> j  iff   i < j
        
    """
    A = {}
    for i in range(5):
        A[i] = []
    
    for i in range(1,5):
        A[0].append(i)
    
    A[1].append(2)
    A[3].append(4)
    
    return A


def IndexToVector( S, index ):
    """
    Given an index, return the configuration in (2S+1)-inary
    """
    length = 5
    N = int( round( (2*S + 1) ) )
    vec = []
    x = int(index)
    while (x != 0):
        vec.append( x % N )
        x /= N
    
    if (length != None):
        if length < len( vec ):
            raise ValueError, "Suggested length is less than the length of the vec number"
        else:
            vec += [0]*(length - len(vec))
    
    vec.reverse()
    
    return vec


def VectorToIndex( S, vec ):
    """
    Given a vector in (2S+1)-inary, return the index
    """
    N = int( round( (2*S + 1) ) )
    return scipy.sum( [ N**i * vec[-(i+1)] for i in range(len(vec)) ] )


def Swap( vec, i, j ):
    """
    Swap elements i and j in the vector vec
    """
    retVec = list(vec)
    retVec[i] = vec[j]
    retVec[j] = vec[i]
    return retVec





##########################################
##                                      ##
##  Simulating various parameter sets   ##
##                                      ##
##########################################


def RunTrials( paramList ):
    """
    Given a list of paramaters of the type
    
        [ [S, Jx, Jz, d], [S2, Jx2, Jz2, d2], ... ]
    
    run the trials, and save them in the BowTieCheck directory
    """
    for paramSet in paramList:
        RunTrial( paramSet )


def RunTrial( paramSet ):
    """
    For just one trial, run and save it
    """
    print '\t' + str(paramSet)
    S, Jx, Jz, d = paramSet
    w, v = BowTieEigenSys( S, Jx, Jz, d )
    
    fileRoot = ParamPath( paramSet )
    
    eigValFile = fileRoot + 'eigVals.pkl'
    eigVecFile = fileRoot + 'eigVecs.pkl'
    
    with open( eigValFile, 'wb' ) as f:
        pickle.dump( w, f )
    with open( eigVecFile, 'wb' ) as f:
        pickle.dump( v, f )


def RunTrials_paramSet( paramSet, printing = False ):
    """
    For many trials, run and save them
    """
    for (w, v, params) in BowTieEigenSys_paramSet( paramSet ):
    
        S, Jx, Jz, d = params
        if printing:
            print '\t' + str( params )
        fileRoot = ParamPath( params )
        
        eigValFile = fileRoot + 'eigVals.pkl'
        eigVecFile = fileRoot + 'eigVecs.pkl'
        
        with open( eigValFile, 'wb' ) as f:
            pickle.dump( w, f )
        with open( eigVecFile, 'wb' ) as f:
            pickle.dump( v, f )



######################################
##                                  ##
##  Loading finished simulations    ##
##                                  ##
######################################

def LoadTrial( paramSet ):
    """
    Given a parameter set [S, Jx, Jz, d], unpickle the eigenvectors and
    eigenvalues (if they exist)
    """
    return LoadEigVals( paramSet ), LoadEigVecs( paramSet )


def LoadEigVals( paramSet ):
    fileRoot = ParamPath( paramSet )
    
    eigValFile = fileRoot + 'eigVals.pkl'
    
    if not os.path.isfile( eigValFile ):
        raise ValueError, "Trial hasn't been recorded"
    else:
        with open( eigValFile, 'rb' ) as f:
            w = pickle.load( f )
        return w

def LoadEigVecs( paramSet ):
    fileRoot = ParamPath( paramSet )
    
    eigVecFile = fileRoot + 'eigVecs.pkl'
    
    if not os.path.isfile( eigVecFile ):
        raise ValueError, "Trial hasn't been recorded"
    else:
        with open( eigVecFile, 'rb' ) as f:
            v = pickle.load( f )
        return v



##########################
##                      ##
##  Other Calculations  ##
##                      ##
##########################

def GroundStateDegeneracy( paramSet, printing = False, threshold = 1e-12 ):
    """
    Given a set of parameters, determine the degeneracy of the ground state
    """
    try:
        w, v = LoadTrial( paramSet )
        groundStateEnergy = w.min()
        indexList = []
        for i in range(len(w)):
            if abs(w[i] - groundStateEnergy) < threshold:
                indexList.append(i)
        if printing:
            print "\tGround State Energy is " + str( groundStateEnergy )
            print "\tThere is/are " + str(len(indexList)) + " of these states"
        return w[indexList], v[:,indexList]
    except ValueError:
        calc = raw_input( 'Not calculated: should we? (y or n)\t' )
        if calc == 'y':
            RunTrial( paramSet )
            return GroundStateDegeneracy( paramSet, threshold )


def GroundStateMeasurement( measurementMatrices, paramSet, printing = False, threshold = 1e-12 ):
    """
    Unpickle the info in paramSet, find the degenerate states, and calculate
    < M > for those states
    """
    w, v = GroundStateDegeneracy( paramSet, printing, threshold )
    
    return w, v, [ [scipy.dot( v[:,i], scipy.dot( measurementMatrix, v[:,i]) ) for i in range(len(w))] for measurementMatrix in measurementMatrices ]
    

def GroundStateSz( S ):
    """
    Calculate < Sz > for the (possible degenerate) ground state energy in
    paramSet
    """
    N = int( round( 2*S + 1 ) )**5
    Sz = []
    for i in range(5):
        Sz.append( scipy.zeros( N ) )
    for i in range(5):
        for index in range(N):
            Sz[i][index] += ( IndexToVector(S, index)[i] - S )
    
    Sz = [scipy.diag(Sz[i]) for i in range(5)]

    return Sz


def GroundStateSz2( S ):
    """
    Calculate < Sz ^ 2 > for the (possible degenerate) ground state energy in
    paramSet
    """
    N = int( round( 2*S + 1 ) )**5
    Sz2 = []
    for i in range(5):
        Sz2.append( scipy.zeros( N ) )
    for i in range(5):
        for index in range(N):
            Sz2[i][index] += ( IndexToVector(S, index)[i] - S )**2.
    
    Sz2 = [scipy.diag(Sz2[i]) for i in range(5)]

    return Sz2


def GroundStateSPlus( S ):
    """
    Calculate < S+_i > for the (possible degenerate) ground state energy in
    paramSet
    """
    N = int( round( 2*S + 1 ) )**5
    SPlus = []
    for i in range(5):
        SPlus.append( scipy.zeros( (N, N) ) )
    for i in range(5):
        for index in range(N):
            config = IndexToVector(S, index)
            sI = config[i]
            if (sI != int( round(2*S) ) ):
                copyConfig = list(config)
                copyConfig[i] += 1
                copyIndex = VectorToIndex( S, copyConfig )
                SPlus[i][copyIndex, index] += scipy.sqrt( (2*S - sI) * (sI + 1) )

    return SPlus


def GroundStateSMinus( S ):
    """
    Calculate < S-_i > for the (possible degenerate) ground state energy in
    paramSet
    """
    N = int( round( 2*S + 1 ) )**5
    SMinus = []
    for i in range(5):
        SMinus.append( scipy.zeros( (N, N) ) )
    for i in range(5):
        for index in range(N):
            config = IndexToVector(S, index)
            sI = config[i]
            if (sI != 0):
                copyConfig = list(config)
                copyConfig[i] -= 1
                copyIndex = VectorToIndex( S, copyConfig )
                SMinus[i][copyIndex, index] += scipy.sqrt( ( sI ) * (2*S - sI + 1) )

    return SMinus


def UpdateGroundStateSCalcs( S, printing = False ):
    """
    Make sure that all recorded values of d are accounted for in the
    'DegenerateSMeasurements.pkl' object
    """
    paramList = ParamList( S )
    SMeasurements = GroundStateSz( S ) + GroundStateSz2( S ) + GroundStateSPlus( S ) + GroundStateSMinus( S )
    with open( 'DegenerateSMeasurements.pkl', 'rb' ) as f:
        masterDic = pickle.load( f )
    for param in paramList:
        if not tuple(param) in masterDic:
            S1, Jx, Jz, d = param
            if printing:
                print '\t' + str( param )
            masterDic[ tuple(param) ] = GroundStateMeasurement( SMeasurements, param )
    with open( 'DegenerateSMeasurements.pkl', 'wb' ) as f:
        pickle.dump( masterDic, f )


def PlotGroundStateExpectations_dSlice( S = 1.0, Jx = 1.0, Jz = 1.0 ):
    """
    Scatter plot of the expectations < Sz > and < Sz^2 > for all the values of
    our recorded trials with S, Jx, and Jz.
    """
    with open('DegenerateSMeasurements.pkl', 'rb') as f:
        masterDic = pickle.load(f)
    
    #   Make the slice.
    for key in masterDic.keys():
        
        if key[:3] != (S, Jx, Jz):
            masterDic.pop( key )
    
    #   Collect a sorted list of the d values
    dList = scipy.array(sorted(masterDic.keys()))[:,-1]
    
    #   Scatter Plot time mutha truckas
    x = []
    E = []
    sZ = [[], [], [], [], []]
    sZ2 = [[], [], [], [], []]
    sPlus = [[], [], [], [], []]
    sMinus = [[], [], [], [], []]

    for d in dList:
        key = (S, Jx, Jz, d)
        wList, vList, sList = masterDic[key]
        
        E.append( wList.min() )	
        
        SzList = sList[:5]
        Sz2List = sList[5:10]
        SPlusList = sList[10:15]
        SMinusList = sList[15:]
        
        l = len(wList)
        
        x += [d]*l
        
        for i in range(5):
            sZ[i]  += SzList[i]
            sZ2[i] += Sz2List[i]
            sPlus[i] += SPlusList[i]
            sMinus[i] += SMinusList[i]
    
    #   E
    pylab.figure(0)
    pylab.plot( dList, E, 'b.' )
    pylab.xlabel('$d$')
    pylab.ylabel('$E$')
    pylab.title('Ground State Energy vs. $d$')
    pylab.legend()
    
    #   dE / dJx
    E = scipy.array( E )
    pylab.figure(-1)
    delE = (E[1:] - E[:-1]) / (dList[1:] - dList[:-1])
    deld = .5 * (dList[1:] + dList[:-1])
    pylab.plot( deld, delE, 'b.' )
    pylab.xlabel('$d$')
    pylab.ylabel('$E$')
    pylab.title('$\partial_{d} E_{GS}$ vs. $d$')
    pylab.legend()
    
    #   Sz
    pylab.figure(1)
    pylab.scatter( x, sZ[0], label = "site 0", color = 'blue' )
    pylab.scatter( x, sZ[1], label = "site 1", color = 'green' )
    pylab.scatter( x, sZ[2], label = "site 2", color = 'red' )
    pylab.scatter( x, sZ[3], label = "site 3", color = 'cyan' )
    pylab.scatter( x, sZ[4], label = "site 4", color = 'magenta' )
    pylab.xlabel('$d$')
    pylab.ylabel('$\left\langle S_i^z \\right\\rangle$')
    pylab.title('$\left\langle \Psi_{GS} \\right| S_i^z \left| \Psi_{GS} \\right\\rangle$ vs. $d$')
    pylab.legend()

    #   Sz2
    pylab.figure(2)
    pylab.scatter( x, sZ2[0], label = "site 0", color = 'blue' )
    pylab.scatter( x, sZ2[1], label = "site 1", color = 'green' )
    pylab.scatter( x, sZ2[2], label = "site 2", color = 'red' )
    pylab.scatter( x, sZ2[3], label = "site 3", color = 'cyan' )
    pylab.scatter( x, sZ2[4], label = "site 4", color = 'magenta' )
    pylab.xlabel('$d$')
    pylab.ylabel('$\left\langle (S_i^z)^2 \\right\\rangle$')
    pylab.title('$\left\langle \Psi_{GS} \\right| (S_i^z)^2 \left| \Psi_{GS} \\right\\rangle$ vs. $d$')
    pylab.legend()
    
    #   SPlus
    pylab.figure(3)
    pylab.scatter( x, sPlus[0], label = "site 0", color = 'blue' )
    pylab.scatter( x, sPlus[1], label = "site 1", color = 'green' )
    pylab.scatter( x, sPlus[2], label = "site 2", color = 'red' )
    pylab.scatter( x, sPlus[3], label = "site 3", color = 'cyan' )
    pylab.scatter( x, sPlus[4], label = "site 4", color = 'magenta' )
    pylab.xlabel('$d$')
    pylab.ylabel('$\left\langle S_i^+ \\right\\rangle$')
    pylab.title('$\left\langle \Psi_{GS} \\right| S_i^+ \left| \Psi_{GS} \\right\\rangle$ vs. $d$')
    pylab.legend()
    
    #   SMinus
    pylab.figure(4)
    pylab.scatter( x, sMinus[0], label = "site 0", color = 'blue' )
    pylab.scatter( x, sMinus[1], label = "site 1", color = 'green' )
    pylab.scatter( x, sMinus[2], label = "site 2", color = 'red' )
    pylab.scatter( x, sMinus[3], label = "site 3", color = 'cyan' )
    pylab.scatter( x, sMinus[4], label = "site 4", color = 'magenta' )
    pylab.xlabel('$d$')
    pylab.ylabel('$\left\langle S_i^- \\right\\rangle$')
    pylab.title('$\left\langle \Psi_{GS} \\right| S_i^- \left| \Psi_{GS} \\right\\rangle$ vs. $d$')
    pylab.legend()
    
    #   Setting the ranges correctly
    dRange = (dList.min(), dList.max())
    for i in range(1,5):
        pylab.figure(i)
        pylab.xlim(dRange)


def PlotGroundStateExpectations_JxSlice( S = 1.0, Jz = 1.0, d = 1.0 ):
    """
    Scatter plot of the expectations < Sz > and < Sz^2 > for all the values of
    our recorded trials with S, Jx, and d
    """
    with open('DegenerateSMeasurements.pkl', 'rb') as f:
        masterDic = pickle.load(f)
    
    #   Make the slice.
    for key in masterDic.keys():
        
        if (key[0] != S) or (key[2:] != (Jz, d)):
            masterDic.pop( key )
    
    #   Collect a sorted list of the Jx values
    JxList = scipy.array(sorted(masterDic.keys()))[:,1]
    
    #   Scatter Plot time mutha truckas
    x = []
    E = []
    sZ = [[], [], [], [], []]
    sZ2 = [[], [], [], [], []]
    sPlus = [[], [], [], [], []]
    sMinus = [[], [], [], [], []]

    for Jx in JxList:
        key = (S, Jx, Jz, d)
        wList, vList, sList = masterDic[key]
        
        E.append( wList.min() )	
        
        SzList = sList[:5]
        Sz2List = sList[5:10]
        SPlusList = sList[10:15]
        SMinusList = sList[15:]
        
        l = len(wList)
        
        x += [Jx]*l
        
        for i in range(5):
            sZ[i]  += SzList[i]
            sZ2[i] += Sz2List[i]
            sPlus[i] += SPlusList[i]
            sMinus[i] += SMinusList[i]
    
    #   E
    pylab.figure(0)
    pylab.plot( JxList, E, 'b.' )
    pylab.xlabel('$J_x$')
    pylab.ylabel('$E$')
    pylab.title('Ground State Energy vs. $J_x$')
    pylab.legend()
    
    #   dE / dJx
    E = scipy.array( E )
    pylab.figure(-1)
    delE = (E[1:] - E[:-1]) / (JxList[1:] - JxList[:-1])
    delJx = .5 * (JxList[1:] + JxList[:-1])
    pylab.plot( delJx, delE, 'b.' )
    pylab.xlabel('$J_x$')
    pylab.ylabel('$E$')
    pylab.title('$\partial_{J_x} E_{GS}$ vs. $J_x$')
    pylab.legend()
    
    #   Sz
    pylab.figure(1)
    pylab.scatter( x, sZ[0], label = "site 0", color = 'blue' )
    pylab.scatter( x, sZ[1], label = "site 1", color = 'green' )
    pylab.scatter( x, sZ[2], label = "site 2", color = 'red' )
    pylab.scatter( x, sZ[3], label = "site 3", color = 'cyan' )
    pylab.scatter( x, sZ[4], label = "site 4", color = 'magenta' )
    pylab.xlabel('$J_x$')
    pylab.ylabel('$\left\langle S_i^z \\right\\rangle$')
    pylab.title('$\left\langle \Psi_{GS} \\right| S_i^z \left| \Psi_{GS} \\right\\rangle$ vs. $J_x$')
    pylab.legend()

    #   Sz2
    pylab.figure(2)
    pylab.scatter( x, sZ2[0], label = "site 0", color = 'blue' )
    pylab.scatter( x, sZ2[1], label = "site 1", color = 'green' )
    pylab.scatter( x, sZ2[2], label = "site 2", color = 'red' )
    pylab.scatter( x, sZ2[3], label = "site 3", color = 'cyan' )
    pylab.scatter( x, sZ2[4], label = "site 4", color = 'magenta' )
    pylab.xlabel('$J_x$')
    pylab.ylabel('$\left\langle (S_i^z)^2 \\right\\rangle$')
    pylab.title('$\left\langle \Psi_{GS} \\right| (S_i^z)^2 \left| \Psi_{GS} \\right\\rangle$ vs. $J_x$')
    pylab.legend()
    
    #   SPlus
    pylab.figure(3)
    pylab.scatter( x, sPlus[0], label = "site 0", color = 'blue' )
    pylab.scatter( x, sPlus[1], label = "site 1", color = 'green' )
    pylab.scatter( x, sPlus[2], label = "site 2", color = 'red' )
    pylab.scatter( x, sPlus[3], label = "site 3", color = 'cyan' )
    pylab.scatter( x, sPlus[4], label = "site 4", color = 'magenta' )
    pylab.xlabel('$J_x$')
    pylab.ylabel('$\left\langle S_i^+ \\right\\rangle$')
    pylab.title('$\left\langle \Psi_{GS} \\right| S_i^+ \left| \Psi_{GS} \\right\\rangle$ vs. $J_x$')
    pylab.legend()
    
    #   SMinus
    pylab.figure(4)
    pylab.scatter( x, sMinus[0], label = "site 0", color = 'blue' )
    pylab.scatter( x, sMinus[1], label = "site 1", color = 'green' )
    pylab.scatter( x, sMinus[2], label = "site 2", color = 'red' )
    pylab.scatter( x, sMinus[3], label = "site 3", color = 'cyan' )
    pylab.scatter( x, sMinus[4], label = "site 4", color = 'magenta' )
    pylab.xlabel('$J_x$')
    pylab.ylabel('$\left\langle S_i^- \\right\\rangle$')
    pylab.title('$\left\langle \Psi_{GS} \\right| S_i^- \left| \Psi_{GS} \\right\\rangle$ vs. $J_x$')
    pylab.legend()
    
    #   Setting the ranges correctly
    JxRange = (JxList.min(), JxList.max())
    for i in range(1,5):
        pylab.figure(i)
        pylab.xlim(JxRange)


def MultiRange( boundSet, stepSet ):
    """
    Given an N x 2 list of range boundaries, and a single step size or an N-long
    list of step sizes, glue together ranges in to cover those bounds in steps
    of that (those) size(s)
    """
    retVal = []
    for i in range( len( boundSet ) ):
        
        bounds = boundSet[ i ]
        if type( stepSet ) != float:
            step = stepSet[ i ]
        else:
            step = stepSet
        
        x = bounds[0]
        xMax = bounds[1] + step
        
        while x < xMax:
            retVal.append( x )
            x += step
    
    return scipy.array( retVal )



##########################
##                      ##
##  File Navigation     ##
##                      ##
##########################


def ParamPath( paramSet ):
    S, Jx, Jz, d = paramSet
    fileRoot = 'Results/S_%(S)03.1f_Jx_%(Jx)03.3f_Jz_%(Jz)03.3f_d_%(d)f_' % {'S':S, 'Jx':Jx, 'Jz':Jz, 'd':d}
    return fileRoot


def ParamList( S ):
    """
    Return every one of the parameter sets we've recorded at this point for spin
    value S
    """
    fileList = os.listdir( 'Results/' )
    fileList = [f for f in fileList if f.startswith( 'S_%(S)03.1f' % {'S':S} ) and f.endswith('_eigVecs.pkl') ]
    fileList = scipy.array([f.split('_') for f in fileList])
    fileList = fileList[:,[1,3,5,7]]
    fileList = [ [ float(s) for s in row ] for row in fileList ]
    fileList = sorted( fileList )
    return fileList
