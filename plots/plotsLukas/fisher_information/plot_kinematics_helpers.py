''' Analysis script for standard plots
'''

# Standard imports and batch mode
import ROOT, os, itertools
ROOT.gROOT.SetBatch(True)
import numpy as np

# RootTools
from RootTools.core.standard             import *
from TTXPheno.Tools.WeightInfo           import WeightInfo

# Import helpers
from plot_helpers                        import *

# Import additional
from array                               import array


def getFisherKinematicsHisto(sample, order=2, var='genZ_pt', plotstring='p_{T}(Z) [GeV]', binning=[10,0,500], selectionString=None, weightString=1., variables=None, parameterList=[]):

    # Polynomial parametrization
    w = WeightInfo(sample.reweight_pkl)
    w.set_order(int(order))

    if variables is None: variables = w.variables

    #remove initial selection string
    sample.setSelectionString('1')
    coeffList = get2DCoeffPlotFromDraw( sample, order, var, binning, selectionString, weightString=weightString )
    detILists  = [ [ np.linalg.det( w.get_fisherInformation_matrix( coeffs, variables, **param['WC'] )[1] ) for coeffs in coeffList ] for param in parameterList ]
    expo = 1. / len(variables)

    hists = []

    for detIList in detILists:

        y_graph = array( 'd', [ abs(detI)**expo for detI in detIList ] )

        hist = ROOT.TH1F( 'histo', 'histo', binning[0], binning[1], binning[2] )
        for i in range(binning[0]):
            hist.SetBinContent(i+1, y_graph[i])

        hists.append(hist)

    return hists#[0]
