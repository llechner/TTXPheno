''' plot helper functions
'''

# TTXPheno imports
from TTXPheno.Tools.WeightInfo     import WeightInfo
from TTXPheno.Tools.WeightInfo     import histo_to_list

# Make a coeff histo from a sample
def getCoeffListFromDraw( sample, order, selectionString, weightString = None ):
    ''' Create list of weights using the Draw function
    '''

    # Polynomial parametrization
    w = WeightInfo(sample.reweight_pkl)
    w.set_order(int(order))

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
    w = WeightInfo(sample.reweight_pkl)
    w.set_order(int(order))

    # 2D Plot, Iteration$ is on x
    histo = sample.get2DHistoFromDraw(
        "Iteration$:%s"%variableString,
        [ len(w.combinations), 0, len(w.combinations) ] + binning,
        selectionString = selectionString,
        weightString = 'p_C*(%s)'%weightString if weightString is not None else 'p_C' )

    return [ histo_to_list(histo.ProjectionX("%i_px"%i, i+1, i+1)) for i in range( histo.GetNbinsY() ) ]

