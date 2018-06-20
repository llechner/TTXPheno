''' plot helper functions
'''

# RootTools
from RootTools.core.standard             import *

# TTXPheno imports
from TTXPheno.Tools.WeightInfo     import WeightInfo
from TTXPheno.Tools.WeightInfo     import histo_to_list

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
def getCoeffPlotFromDraw( sample, order, variableString, binning, selectionString, weightString = None ):
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

    return [ histo_to_list( histo.ProjectionX("%i_px"%i, i+1, i+1) ) for i in range( histo.GetNbinsY() ) ]


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
        coeffs.append( [ reader.event.p_C[i] * weightFunction( reader.event, sample ) if weightFunction is not None else reader.event.p_C[i] for i in range(reader.event.np) ] )

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

