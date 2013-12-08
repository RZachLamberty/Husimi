################################################################################
################################################################################
##                                                                            ##
##  CorrelationScaling.py                                                     ##
##                                                                            ##
##      This is a collection of routines designed to measure the decay        ##
##      scaling of correlation functions measured in the infinite Husimi      ##
##      cactus DMRG code.  We claim LRO is indicated by decays in these       ##
##      correlations of the form 2^-R, so we will fit our correlations to     ##
##      the functional form A[0] * A[1]**-R                                   ##
##                                                                            ##
################################################################################
################################################################################

import scipy
import scipy.optimize as so
import pylab
import cPickle as pickle
from Plotter import R, SZ, SZ2, SIZSCZ, SIPSCM, SIMSCP, SIDOTSC, LABEL_DIC, BEST_GEN



#########################
#                       #
#   Module constants    #
#                       #
#########################

MEASUREMENT_INFO_FILE = 'TrialArrayDic.pkl'
SCALING_INFO_FILE = 'ScalingConstantsDic.pkl'

CORRELATION_INDICES = [ SIZSCZ, SIPSCM, SIMSCP, SIDOTSC ]

for key in LABEL_DIC:
    LABEL_DIC[key] = LABEL_DIC[key].replace('$', '')


#############################
#                           #
#   Performing the fits     #
#                           #
#############################

def AllScalingInfo( fitStyle = 'powerLaw', userRange = False ):
    """
    Cycle through the valid correlation functions and see if there is any
    reasonable fit to the form A[0] * A[1]^(-R)
    """
    
    with open( MEASUREMENT_INFO_FILE, 'rb' ) as f:
        x = pickle.load( f )
    
    with open( SCALING_INFO_FILE, 'rb' ) as f:
        scalingDic = pickle.load( f )
    
    for trial in x:
        
        print trial
        if trial not in scalingDic:
            
            scalingDic[trial] = {}
            
            for gen in x[trial]:
                
                #   We don't start fitting until R = 2, so let's at least have some
                #   points to fit!
                if gen >=4:
                    scalingDic[trial][gen] = {}
                    
                    if userRange:
                        fig = pylab.figure(-2)
                        p = fig.add_subplot(111)
                        p.plot( x[trial][gen][0][:,R], abs(x[trial][gen][0][:,max(CORRELATION_INDICES)]), linestyle = '', marker = 'o', ms = 10 )
                        p.set_yscale( 'log' )
                        lv = raw_input('Low Val?    ')
                        if lv != '':
                            lv = int(eval(lv))
                        else:
                            lv = 2
                        hv = raw_input('High Val?   ')
                        if hv != '':
                            hv = int(eval(hv)) + 1
                        else:
                            hv = len(x[trial][gen][0])
                        fig.clf()
                    else:
                        lv, hv = 2, len(x[trial][gen][0])
                    
                    for sz in x[trial][gen]:
                        
                        A = x[trial][gen][sz].copy()
                        
                        scalingDic[trial][gen][sz] = FitArray( A, fitStyle, (lv, hv) )
    
    
    with open( SCALING_INFO_FILE, 'wb' ) as f:
        pickle.dump( scalingDic, f )
    

def FitArray( A, fitStyle = 'powerLaw', userRange = (0,1) ):
    """
    Given an array of an arbitrary size, fit it to the functional form described
    by function FitFunc and CostFunc.  The information gained in that fitting is
    returned in the form
    
        { 'Correlation label' : ( arrayA, chiSquared ), ... }
    """
    retDic = {}
    
    #   First, get everything in proper order:
    if A[0,0] != 0.0:
        A[:,R] = A[:,R].max() - A[:,R]
        orderIndices = A[:,R].argsort()
        A = A[orderIndices]
    
    #   Now perform fits, and return the fitting information in a dictionary
    #if userRange:
    #    fig = pylab.figure(-2)
    #    p = fig.add_subplot(111)
    #    p.plot( A[:,R], abs(A[:,max(CORRELATION_INDICES)]), linestyle = '', marker = 'o', ms = 10 )
    #    p.set_yscale( 'log' )
    #    lv = int( eval( raw_input('Low Val?    ') ) )
    #    hv = int( eval( raw_input('High Val?   ') ) ) + 1
    #    fig.clf()
    #else:
    #    lv = 2
    #    hv = len(A)
    lv, hv = userRange
        
    x = A[lv:hv,R].copy()
    for index in CORRELATION_INDICES:
        y = abs( A[lv:hv, index].copy() )
        
        if fitStyle == 'powerLaw':
            aGuess = scipy.array([1.,1.])
            aTup = so.fmin( costFunc, aGuess, (x, y), full_output = True, disp = False )
        elif fitStyle == 'log':
            aGuess = scipy.array([0.,1.])
            aTup = so.fmin( costFuncLine, aGuess, (x, scipy.log(y)), full_output = True, disp = False )
        else:
            raise ValueError, "Not a valid type of fit yet.  Try powerLaw or log"
        
        aOptimal, chiSquared, iterations, funcalls, warnflag = aTup
        
        if warnflag == 0:
            if fitStyle == 'powerLaw':
                retDic[ LABEL_DIC[ index ] ] = aOptimal, chiSquared
            elif fitStyle == 'log':
                retDic[ LABEL_DIC[ index ] ] = scipy.exp(aOptimal), chiSquared
    
    return retDic


def FitFunc( a, x ):
    return a[0] * a[1]**scipy.array(-x)


def costFunc( a, x, y ):
    return scipy.sum( ( FitFunc( a, x ) - y)**2. )


def FitFuncLine( a, x ):
    return a[0] - a[1] * x


def costFuncLine( a, x, y ):
    return scipy.sum( ( FitFuncLine( a, x ) - y)**2. )




#############################
#                           #
#   Comparing the fits      #
#                           #
#############################

def CompareFits( trialList = [] ):
    """
    Cycle through the trials; pick out the best generation, and see what fits
    we came up with
    """
    
    with open( MEASUREMENT_INFO_FILE, 'rb' ) as f:
        dataDic = pickle.load( f )
    
    with open( SCALING_INFO_FILE, 'rb' ) as f:
        fitDic = pickle.load( f )
    
    colorList = [ 'blue', 'green', 'red', 'cyan', 'magenta', 'gold', 'black' ]
    lineStyleList = [ '-', '--', '-.', '..' ]
    markerList = [ 'o', '^', '+', 'D', 'v', 'x', 's' ]
    
    if trialList == []:
        trialList = fitDic.keys()
    
    for trial in trialList:
        bestGen = BEST_GEN[ trial ]
        
        #   In case bestGen wasn't in, cycle downward until it is
        while (bestGen > 4) and (bestGen not in fitDic[trial]):
            bestGen -= 1
        
        thisFitDic = fitDic[ trial ][ bestGen ]
        thisDataDic = dataDic[ trial ][ bestGen ]
        szList = sorted( thisFitDic.keys() )
        
        for i in range( len( CORRELATION_INDICES ) ):
            
            corIndex = CORRELATION_INDICES[i]
            
            fig = pylab.figure( i )
            fig.clf()
            p = fig.add_subplot( 111 )
        
            for j in range( len( szList ) ):
                
                sz = szList[j]
                
                if LABEL_DIC[corIndex] in thisFitDic[sz]:
                    
                    #   Probably don't need individual sub-plots...
                    #p = fig.add_subplot( len( szList ), 1, (j + 1) )
                    
                    c = colorList[j]
                    
                    #   Data -- plotted with solid markers
                    A = thisDataDic[sz].copy()
                    if A[0,0] != 0.0:
                        A[:,R] = A[:,R].max() - A[:,R]
                        o = A[:,R].argsort()
                        A = A[o]
                    x = A[:,R]
                    y = A[:,corIndex]
                    p.plot( x, abs(y), linestyle = '', marker = 'o', mec = 'k', mfc = c, mew = 2, ms = 9, label = '$S_z = ' + str(sz) + '$' )
                    
                    #   Fit -- dashed lines
                    aTup = thisFitDic[sz][ LABEL_DIC[ corIndex ] ]
                    a, cost = aTup
                    p.plot( x, FitFunc( a, x ), linestyle = '--', linewidth = 3, color = c, label = '$a_1 = ' + str(a[1]) + '$' )
            
            yl = p.get_ylim()
            p.set_yscale('log')
            p.set_ylim( yl )
            p.set_xlabel( '$R$' )
            p.set_ylabel( '$' + LABEL_DIC[ corIndex ] + '$' )
            p.set_title( '$' + LABEL_DIC[ corIndex ] + '$' )
            p.legend(loc = 'lower left')
        
        raw_input( str( trial ) )
    
    pylab.close( 'all' )


def TrendsInConstants( trialList = [] ):
    """
    Cycle through the trials; Plot the values obtained for constants as a
    function of the generation
    """
    
    with open( SCALING_INFO_FILE, 'rb' ) as f:
        fitDic = pickle.load( f )
    
    colorList = [ 'blue', 'green', 'red', 'cyan', 'magenta', 'gold', 'black' ]
    lineStyleList = [ '-', '--', '-.', '..' ]
    markerList = [ 'o', '^', '+', 'D', 'v', 'x', 's' ]
    
    if trialList == []:
        trialList = fitDic.keys()
    
    for trial in trialList:
    
        try:
            genList = [ gen for gen in fitDic[trial] if ( (len(fitDic[trial][gen][0]) != 0) and ('S_i \\cdot S_c' in fitDic[trial][gen][0]) ) ]
            aValList = [ fitDic[trial][gen][0]['S_i \\cdot S_c'][0][1] for gen in genList ]
            
            fig = pylab.figure(-2)
            p = fig.add_subplot(111)
            p.plot( genList, aValList, linestyle = '', marker = 'o', ms = 10 )
            p.vlines( BEST_GEN[trial], 2., 3. )
            p.set_ylim((2., 2.5))
            
            raw_input( str( (trial, aValList[genList.index(BEST_GEN[trial])] ) ) )
            pylab.clf()
        except:
            raw_input( str( trial ) + ", something went wrong..." )
    
    pylab.close( 'all' )
