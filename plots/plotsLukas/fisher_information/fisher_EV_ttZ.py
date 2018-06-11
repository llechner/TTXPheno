#!/usr/bin/env python
''' Analysis script for standard plots
'''

# Standard imports and batch mode
import ROOT, os, itertools
ROOT.gROOT.SetBatch(True)
from math                                import sqrt, cos, sin, pi, isnan, sinh, cosh
import copy
import imp
import numpy as np

# RootTools
from RootTools.core.standard             import *

# TTXPheno
from TTXPheno.Tools.user                 import plot_directory
from TTXPheno.Tools.helpers              import deltaPhi, getCollection, deltaR 
from TTXPheno.Tools.WeightInfo           import WeightInfo
from TTXPheno.Tools.cutInterpreter       import cutInterpreter

# Import samples
from TTXPheno.samples.benchmarks         import *

# Import helpers
from plot_helpers                        import *

# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--plot_directory',     action='store',      default='gen')
argParser.add_argument('--sample',             action='store',      default='fwlite_ttZ_ll_LO_order3_8weights')
argParser.add_argument('--order',              action='store',      default=2)
argParser.add_argument('--selection',          action='store',      default='lepSel3-onZ-njet3p-nbjet1p-Zpt0', help="Specify cut.")
argParser.add_argument('--small',              action='store_true', help='Run only on a small subset of the data?')

argParser.add_argument('--parameters',         action='store',      default = ['ctW', '3', 'ctWI', '3', 'ctZ', '3', 'ctZI', '3'], type=str, nargs='+', help = "argument parameters")
argParser.add_argument('--luminosity',         action='store',      default=150)

args = argParser.parse_args()

# Import samples
sample_file = "$CMSSW_BASE/python/TTXPheno/samples/benchmarks.py"
samples = imp.load_source( "samples", os.path.expandvars( sample_file ) )
sample = getattr( samples, args.sample )

# Scale the plots with number of events used (implemented in ref_lumiweight1fb)
event_factor = 1.
if args.small:
    sample.reduceFiles( to = 1 )
    event_factor = sample.nEvents / 5000.

# Polynomial parametrization
w = WeightInfo(sample.reweight_pkl)
w.set_order(int(args.order))

selection_string = cutInterpreter.cutString( args.selection )
weightString = 10000 #'ref_lumiweight1fb*%s*%s' %( str(args.luminosity), str(event_factor) )

variables = [
#             'cpQM',
#             'cpt',
             'ctZ',
#             'ctZI'
]

plot_variables = [
                  { 'name':'Z_pt', 'binning':[ 20, 0, 500 ] },
                  { 'name':'Z_mass', 'binning':[ 20, 70, 110 ] },
                  { 'name':'Z_phi', 'binning':[ 20, -np.pi, np.pi ] },
                  { 'name':'Z_eta', 'binning':[ 20, -3, 3 ] },
                  { 'name':'Z_cosThetaStar', 'binning':[ 20, -1.2, 1.2 ] },
]

coeff_tot = getCoeffListFromDraw( sample, args.order, None, weightString=weightString )
coeff_sel = getCoeffListFromDraw( sample, args.order, selection_string, weightString=weightString )
coeff_var = [ getCoeffPlotFromDraw( sample, args.order, var['name'], var['binning'], selection_string, weightString=weightString ) for var in plot_variables ]

detI_tot = np.linalg.det( w.get_fisherInformation_matrix( coeff_tot, variables )[1] )
detI_sel = np.linalg.det( w.get_fisherInformation_matrix( coeff_sel, variables )[1] )

detI_var = [ np.linalg.det( w.get_total_fisherInformation_matrix( coeff, variables )[1] ) for coeff in coeff_var ]

det = np.array([detI_tot, detI_sel] + detI_var)# / detI_tot

print( len(det) )
print( np.range( 1, len(det)+1 ) )
print( det )

c1 = ROOT.TCanvas()
gr = ROOT.TGraph( len(det), range( 1, len(det)+1 ), det )
#gr.SetLineColor( 2 )
#gr.SetLineWidth( 4 )
#gr.SetMarkerColor( 4 )
#gr.SetMarkerStyle( 21 )
#gr.SetTitle( 'a simple graph' )
#gr.GetXaxis().SetTitle( 'X title' )
#gr.GetYaxis().SetTitle( 'Y title' )
gr.Draw( 'AB' )
 
 # TCanvas.Update() draws the frame, after which one can change it
c1.Update()
#c1.GetFrame().SetFillColor( 21 )
#c1.GetFrame().SetBorderSize( 12 )
#c1.Modified()
#c1.Update()


plot_directory_ = os.path.join(\
    plot_directory,
    args.plot_directory, 
    sample.name, 
    'fisher_information', 
    args.selection)

c1.Print(os.path.join(plot_directory_, 'fisher_information.png'))

