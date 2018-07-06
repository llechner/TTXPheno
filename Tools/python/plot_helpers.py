''' plot helper functions
'''

# RootTools
from RootTools.core.standard             import *

# TTXPheno imports
from TTXPheno.Tools.WeightInfo     import WeightInfo
from TTXPheno.Tools.WeightInfo     import histo_to_list

import numpy as np
from itertools import chain

# Make a coeff histo from a sample
def getCoeffListFromDraw( sample, order, selectionString, weightString = None ):
    ''' Create list of weights using the Draw function
    '''

    # Polynomial parametrization
    w = WeightInfo( sample.reweight_pkl )
    w.set_order( int(order) )

    # Draw 
    histo = sample.get1DHistoFromDraw(
        "Iteration$",
        [ len(w.combinations), 0, len(w.combinations) ],
        selectionString = selectionString,
        weightString = 'p_C*(%s)'%weightString if weightString is not None else 'p_C' )

    return histo_to_list( histo )

# Make a coeff histo from a sample
def getCoeffPlotFromDraw( sample, order, variableString, binning, selectionString, weightString = None, nEventsThresh = 0 ):
    ''' Create list of weights using the Draw function
    '''

    # Polynomial parametrization
    w = WeightInfo( sample.reweight_pkl )
    w.set_order( int(order) )

    histo = sample.get2DHistoFromDraw(
        "%s:Iteration$"%variableString,
        [ len(w.combinations), 0, len(w.combinations) ] + binning,
        selectionString = selectionString,
        weightString = 'p_C*(%s)' %(weightString) if weightString is not None else 'p_C' )

    return [ histo_to_list( histo.ProjectionX("%i_px"%i, i+1, i+1) ) for i in range( histo.GetNbinsY() ) if histo.ProjectionX("%i_px"%i, i+1, i+1).GetEntries() >= int(nEventsThresh) ]


# Make a coeff histo from a sample
def get2DCoeffPlotFromDraw( sample, order, variableString, binning, selectionString, weightString = None, nEventsThresh = 0 ):
    ''' Create list of weights using the Draw function
    '''

    # Polynomial parametrization
    w = WeightInfo( sample.reweight_pkl )
    w.set_order( int(order) )

    histo = sample.get3DHistoFromDraw(
        "%s:Iteration$"%variableString,
        [ len(w.combinations), 0, len(w.combinations) ] + binning,
        selectionString = selectionString,
        weightString = 'p_C*(%s)' %(weightString) if weightString is not None else 'p_C' )

    return  [ histo_to_list( histo.ProjectionX("%i_%i_px"%(i,j), i+1, i+1, j+1, j+1) ) for i in range( histo.GetNbinsY() ) for j in range( histo.GetNbinsZ() ) if histo.ProjectionX("%i_%i_px"%(i,j), i+1, i+1, j+1, j+1).GetEntries() > int(nEventsThresh) ]


# Make a coeff histo from a sample
def get3DCoeffPlotFromDraw( sample, order, variableString, binning, selectionString, weightString = None, nEventsThresh = 0 ):
    ''' Create list of weights using the Draw function
    '''

    if len(binning) != 9: raise
    if len(variableString.split(':')) != 3: raise

    variableString2D = ':'.join( variableString.split(':')[1:] )
    variableString3D = variableString.split(':')[0]

    bounds = np.linspace( start=binning[7], stop=binning[8], num=binning[6]+1 )
    coeffList3D = []
    for i, bound in enumerate(bounds[1:]):
        sample.setSelectionString('(1)')
        coeffList3D.append( get2DCoeffPlotFromDraw( sample, order, variableString2D, binning[:6], selectionString + '&&%s>=%f&&%s<%f'%( variableString3D, bounds[i], variableString3D, bound ), weightString, nEventsThresh ) )

    coeffList = []
    for coeffs in coeffList3D:
        if len( coeffs ) == 0: continue
        for coeff in coeffs:
            if len( coeff ) == 0: continue
            coeffList.append( coeff )

    return coeffList


def getCoeffListFromEvents( sample, selectionString = None, weightFunction = None ):
    ''' Create list of weights for each event
    '''

    sample.setSelectionString( selectionString ) 

    variables = map( TreeVariable.fromString, [ "np/I", "ref_lumiweight1fb/F", "lumiweight1fb/F" ] )
    variables.append( VectorTreeVariable.fromString('p[C/F]', nMax=1000) )

    reader = sample.treeReader( variables = variables )
    reader.start()

    coeffs = []
    while reader.run():
        coeffs.append( [ reader.event.p_C[i] * (weightFunction( reader.event, sample )) if weightFunction is not None else reader.event.p_C[i] for i in range(reader.event.np) ] )

    return coeffs


def replace_selectionstrings( selectionList ):
    ''' replace selection string elements with string shown in plots
    '''

    replaced_selectionstrings = []

    for item in selectionList:

        # lepton selection
        if 'lepSel' in item:
            replaced_selectionstrings.append(item[-1] + ' lep' if item[-1]=='1' else item[-1] + ' leps')

        # number of jets/bjets/leps selection
        elif item[0] == 'n':
            num = ''.join(c for c in item if c.isdigit())
            particle = item.split('n')[1].split(num)[0]
            replaced_selectionstrings.append( ' '.join([num, particle] if num=='1' else [num, particle + 's']) )

        # pt selection
        elif 'pt' in item:
            particle = item.split('pt')[0]
            val = item.split('pt')[1]

            if 'to' in val:
                replaced_selectionstrings.append( pre + ' #leq p_{T}(' + particle + ') < ' + post + ' GeV')
            else:
               	replaced_selectionstrings.append( 'p_{T}(' + particle + ') > ' + val + ' GeV')

        # anything
        else:
            replaced_selectionstrings.append( item )

    return replaced_selectionstrings

