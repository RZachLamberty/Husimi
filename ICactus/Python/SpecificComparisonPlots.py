################################################################################
################################################################################
##                                                                            ##
##  SpecificComparisonPlots.py                                                ##
##                                                                            ##
##      This is a collection of routines designed to plot specific trial      ##
##      measurements for direct comparison to those found in other studies.   ##
##                                                                            ##
################################################################################
################################################################################

import scipy
import cPickle as pickle
import copy
import pylab
from mpl_toolkits.mplot3d import Axes3D

import Parser
import Plotter as P
from CorrelationScaling import MEASUREMENT_INFO_FILE, SCALING_INFO_FILE
pp = P.Plotter()


################################################################################
#                                                                              #
#   All expectation value plots, broken into regions                           #
#                                                                              #
################################################################################




#############################
#                           #
#   Isakov and Kim plots    #
#                           #
#############################

def DisorderedQuantumParamagnet(twoToRNormalization = False, threeD = False, SzSubset = 'All'):
    """
    Large positive d, small negative Jx
    
    Trials measured in this regime:
        
        (1.0, -0.5, 1.0, 10.0)
        
    """
    
    pp.AllSliceMeasurements( 'Jx', (1.0, 1.0, 10.0, 'best'), TwoToRNormalization = twoToRNormalization, ThreeD = threeD, varRange = (-0.5, 0.0), szSubset = SzSubset )


def XYFerromagnet(twoToRNormalization = False, threeD = False, SzSubset = 'All'):
    """
    Large negative J_x
    
    Trials measured in this regime:
        
        (1, -10.0, 1, 0),
        (1, -0.95, 1, 0),
        (1, -0.90, 1, 0),
        (1, -0.85, 1, 0),
        (1, -0.80, 1, 0),
        (1, -0.75, 1, 0),
        (1, -0.70, 1, 0),
        (1, -0.65, 1, 0),
        (1, -0.60, 1, 0),
        (1, -0.55, 1, 0),
        (1, -0.50, 1, 0),
        (1, -0.45, 1, 0),
        (1, -0.40, 1, 0),
        (1, -0.35, 1, 0),
        (1, -0.30, 1, 0),
        (1, -0.25, 1, 0),
        (1, -0.20, 1, 0),
        (1, -0.15, 1, 0)
    """
    
    if threeD:
        v = (-0.95, -0.11)
    else:
        v = (-10, -0.11)
    
    pp.AllSliceMeasurements( 'Jx', (1.0, 1.0, 0.0, 'best'), TwoToRNormalization = twoToRNormalization, ThreeD = threeD, varRange = v, szSubset = SzSubset )


def Specific_XYFerromagnet(SzSubset = 'All'):
    """
    (1, -10.0, 1, 0),
    """
    
    fixedParams = (1.0, 1.0, 0.0, 'best')
    plotDic, keyList, varList = pp.Jx_Slice( fixedParams )
    
    #   Shed the excess trials (all of them)
    try:
        trial = (1,-10,1,0)
        g = P.BEST_GEN[trial]
        plotDic = plotDic[trial + (g,)]
    except KeyError:
        raise ValueError, "Must have changed the best generation for this trial!"
    
    A = plotDic[0].copy()
    x = A[:,P.R]
    siz = A[:,P.SIZSCZ]
    sPM = A[:,P.SIPSCM]
    sMP = A[:,P.SIMSCP]
    
    fig = pylab.figure()
    p = fig.add_subplot(111)
    p1 = p.plot( x, abs( siz ), c = 'blue', lw = 2, marker = 'o', ms = 10, mew = 2, mec = 'b', mfc = 'w', label = '$\left\langle S_0^z S_R^z \\right \\rangle$' )
    p2 = p.plot( x, sPM, c = 'r', lw = 2, marker = 's', ms = 10, mew = 2, mec = 'r', mfc = 'w', label = '$\left\langle S_0^+ S_R^- \\right \\rangle$' )
    p3 = p.plot( x, sMP, c = 'g', lw = 2, marker = 'd', ms = 10, mew = 2, mec = 'g', mfc = 'w', label = '$\left\langle S_0^- S_R^+ \\right \\rangle$' )
    
    leg = p.legend()
    ltext = leg.get_texts()
    pylab.setp(ltext, fontsize = 26)
    
    p.set_xlabel('$R$', size = 26)
    
    pylab.show()
    

def MysteryRegion(twoToRNormalization = False, threeD = False, SzSubset = 'All'):
    """
    Small negative J_x and very small (or strictly zero) d
    
    Trials measured in this regime:
    
        (1, -0.05, 1, 0.5)
    """
    pp.AllSliceMeasurements( 'Jx', (1.0, 1.0, 0.5, 'best'), TwoToRNormalization = twoToRNormalization, ThreeD = threeD, varRange = (-.5, 0.0), szSubset = SzSubset )
    

def TriangleVBS(twoToRNormalization = False, threeD = False, SzSubset = 'All'):
    """
    Small negative J_x and very small (or strictly zero) d
    
    Trials measured in this regime:
    
        (1, -0.10, 1, 0),
        (1, -0.09, 1, 0),
        (1, -0.05, 1, 0),
        (1, -0.01, 1, 0)
    """
    
    pp.AllSliceMeasurements( 'Jx', (1.0, 1.0, 0.0, 'best'), TwoToRNormalization = twoToRNormalization, ThreeD = threeD, varRange = (-0.1, 0.0), szSubset = SzSubset )




#############################
#                           #
#   Other plots             #
#                           #
#############################

def IsotropicHeisenberg(twoToRNormalization = False, threeD = False, SzSubset = 'All'):
    """
    J_x ~ J_z and d = 0
    """
    
    pp.AllSliceMeasurements( 'Jx', (1.0, 1.0, 0.0, 'best'), TwoToRNormalization = twoToRNormalization, ThreeD = threeD, varRange = (0.5, 1.5), szSubset = SzSubset )



def HeisenbergScalingCorrelation( startVal = 0.5, endVal = 1.5 ):

    #   You may need to flip around some of the hashed values in Plotter.py to
    #   do this up right.
    
    trialList = Parser.KeysInRange( 1, [startVal, endVal], 1, 0 )
    
    with open( SCALING_INFO_FILE, 'rb' ) as f:
        scalingDic = pickle.load( f )
    
    JxList = []
    y = []
        
    for trial in trialList:
        
        bestGen = P.BEST_GEN[ trial ]
        
        #   In case bestGen wasn't in, find the largest dictionary which we can
        #   compare
        #while (bestGen > 4) and ( (bestGen not in scalingDic[trial]) or (len(scalingDic[trial][bestGen][0]) == 0) or ('S_i \\cdot S_c' not in scalingDic[trial][bestGen][0]) ):
        #    bestGen -= 1
        
        if bestGen == 4:
            pass
        else:
            thisVal = scalingDic[ trial ][ bestGen ][ 0 ][ 'S_i \\cdot S_c' ]
            JxList.append( trial[1] )
            y.append( thisVal[0][1] )
    
    #   For Now...
    print "popping..."
    JxList.pop(-2)
    JxList.pop(-2)
    y.pop(-2)
    y.pop(-2)
    
    return JxList, y
    
    #plot( JxList, y, linestyle = '--', marker = 'o', markersize = 7, mfc = 'b', mew = 3, mec = 'k' )


def LargeSpinHeisenbergScalingCorrelation( startVal = 0.9, endVal = 1.1 ):
    
    trialList = Parser.KeysInRange( [1.0,4.0], [startVal, endVal], 1, 0 )
    i = 0
    while i < len(trialList):
        s, jx, jz, d = trialList[i]
        if not ((jx == startVal) or (jx == endVal) or (jx == 1.0)):
            t = trialList.pop(i)
        else:
            i += 1
    
    
    with open( SCALING_INFO_FILE, 'rb' ) as f:
        scalingDic = pickle.load( f )
    
    JxListDic = {}
    yDic = {}
        
    for trial in trialList:
        
        bestGen = P.BEST_GEN[ trial ]
        
        #   In case bestGen wasn't in, find the largest dictionary which we can
        #   compare
        #while (bestGen > 4) and ( (bestGen not in scalingDic[trial]) or (len(scalingDic[trial][bestGen][0]) == 0) or ('S_i \\cdot S_c' not in scalingDic[trial][bestGen][0]) ):
        #    bestGen -= 1
        
        if bestGen <= 4:
            pass
        else:
            thisVal = scalingDic[ trial ][ bestGen ][ 0 ][ 'S_i \\cdot S_c' ]
            if trial[0] not in JxListDic:
                JxListDic[trial[0]] = []
            JxListDic[trial[0]].append( trial[1] )
            if trial[0] not in yDic:
                yDic[trial[0]] = []
            yDic[trial[0]].append( thisVal[0][1] )
    
    return JxListDic, yDic
    
    #
    c = {1:'b', 2:'r', 3:'g', 4:'c'}
    for S in JxListDic.keys():
        JxList = JxListDic[S]
        yList = yDic[S]
        pylab.plot( JxList, yList, linestyle = '--', color = c[S], marker = 'o', markersize = 7, mfc = c[S], mew = 3, mec = 'k' )
