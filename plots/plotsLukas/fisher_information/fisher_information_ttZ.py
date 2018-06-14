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
argParser.add_argument('--sample',             action='store',      default='fwlite_ttZ_ll_LO_order2_15weights_ref')
argParser.add_argument('--order',              action='store',      default=2)
argParser.add_argument('--small',              action='store_true', help='Run only on a small subset of the data?')
argParser.add_argument('--selection',          action='store',      default='lepSel3-onZ-njet3p-nbjet1p-Zpt0', help="Specify cut.")
argParser.add_argument('--variables',          action='store',      default = ['cpt', 'cpQM', 'ctZ', 'ctZI'], type=str, nargs='+', help = "argument parameters")
argParser.add_argument('--parameters',         action='store',      default = ['ctW', '3', 'ctWI', '3', 'ctZ', '3', 'ctZI', '3'], type=str, nargs='+', help = "argument parameters")
argParser.add_argument('--luminosity',         action='store',      default=150)

args = argParser.parse_args()

# Import samples
sample_file = "$CMSSW_BASE/python/TTXPheno/samples/benchmarks.py"
samples = imp.load_source( "samples", os.path.expandvars( sample_file ) )
sample = getattr( samples, args.sample )

# Scale the plots with number of events used (implemented in ref_lumiweight1fb)
event_factor = 1.
fisher_directory = 'fisher_information'
if args.small:
    sample.reduceFiles( to = 1 )
    event_factor = sample.nEvents / float(sample.chain.GetEntries())
    fisher_directory += '_small'

# Polynomial parametrization
w = WeightInfo(sample.reweight_pkl)
w.set_order(int(args.order))

# Format input parameters to dict
WC_string = 'SM'
WC = {}
if len(args.parameters) != 0:
    coeffs = args.parameters[::2]
    str_vals = args.parameters[1::2]
    vals   = list( map( float, str_vals ) )
    WC = { coeff:vals[i] for i, coeff in enumerate(coeffs) }
    WC_string = '_'.join( args.parameters )


def get_reweight_function():
    ''' return a weight function
    '''

    def reweight( event, sample ):
        return event.ref_lumiweight1fb * args.luminosity * event_factor

    return reweight


selection_string = cutInterpreter.cutString( args.selection )

# Make sure that weightString contains the same as weightFunction!!!
weightString = 'ref_lumiweight1fb*%s*%s' %( str(args.luminosity), str(event_factor) )
weightFunction = get_reweight_function()

# split selection string step by step
selectionStrings = []
selectionElements = args.selection.split('-')
for i, item in enumerate(selectionElements[:-1]):
    selectionStrings.append( {'plotstring':' + '.join(replace_selectionstrings(selectionElements[:i+1])) + ' selection', 'selection':'-'.join(selectionElements[:i+1])} )
selectionStrings.append( {'plotstring':'full pre-selection', 'selection':args.selection} )

# Additional plot variables
plot_variables =   [
                     { 'plotstring':'p_{T}(Z)',            'var':'Z_pt',           'binning':[ 20000, 0, 500 ] },
                     { 'plotstring':'p_{T}(Z)',            'var':'Z_pt',           'binning':[ 2000, 0, 500 ] },
                     { 'plotstring':'p_{T}(Z)',            'var':'Z_pt',           'binning':[ 200, 0, 500 ] },
                     { 'plotstring':'p_{T}(Z)',            'var':'Z_pt',           'binning':[ 20, 0, 1000 ] },
                     { 'plotstring':'cos(#theta*)',        'var':'Z_cosThetaStar', 'binning':[ 20, -1.2, 1.2 ] },
                     { 'plotstring':'m_{ll}',              'var':'Z_mass',         'binning':[ 20, 70, 110 ] },
                     { 'plotstring':'#phi(Z)',             'var':'Z_phi',          'binning':[ 20, -np.pi, np.pi ] },
                     { 'plotstring':'#eta(Z)',             'var':'Z_eta',          'binning':[ 20, -3, 3 ] },
                     { 'plotstring':'p_{T}(E_{T}^{miss})', 'var':'GenMet_pt',      'binning':[ 20, 0, 500 ] },
                     { 'plotstring':'#phi(E_{T}^{miss})',  'var':'GenMet_phi',     'binning':[ 20, -np.pi, np.pi ] },
]

# Calculate coefficients for binned distribution
# Calculate determinant using the 'variables' submatrix of FI
for var in plot_variables:
    var['coeff']       = getCoeffPlotFromDraw( sample, args.order, var['var'], var['binning'], selection_string, weightString=weightString )
    var['detI']        = np.linalg.det( w.get_total_fisherInformation_matrix( var['coeff'], args.variables, **WC )[1] )
    # add bin information to plot labels
    var['plotstring'] += ' (' + str(var['binning'][0]) + ' bins)'

# sort by Fisher Information
#plot_variables = sorted(plot_variables, key=lambda k: -k['detI'])

# Calculate coefficients unbinned (event loop)
# Calculate determinant using the 'variables' submatrix of FI
for selection in selectionStrings:
    selection['coeff'] = getCoeffListFromEvents( sample, selectionString = cutInterpreter.cutString(selection['selection']), weightFunction = weightFunction )
    selection['detI']  = np.linalg.det( w.get_total_fisherInformation_matrix( selection['coeff'], args.variables, **WC )[1] )

# sort by Fisher Information
#selectionStrings = sorted(selectionStrings, key=lambda k: -k['detI'])

# Full Fisher information
full          = { 'plotstring':'full'}
full['coeff'] = getCoeffListFromEvents( sample, selectionString = None, weightFunction = weightFunction )
full['detI']  = np.linalg.det( w.get_total_fisherInformation_matrix( full['coeff'], args.variables, **WC )[1] )


# Data for plotting
data = [full] + selectionStrings + plot_variables
expo = 1. / len(args.variables)
detI = [ abs( item['detI'] / full['detI'] )**expo if full['detI'] != 0 else 0 for item in data ]
n = len(detI)

print(detI)

# Prepare data for plotting
from array import array
x = array('d')
y = array('d')

for i in range(n):
    x.append(i+1)
    y.append( detI[i] )

# Plotting
c1 = ROOT.TCanvas()
c1.SetBottomMargin(0.05)

# Plot
gr = ROOT.TGraph( n, x, y )
gr.GetXaxis().SetTickLength(0.)
gr.SetFillColor( 40 )
gr.SetLineWidth( 4 )
gr.GetXaxis().SetLabelSize(0)
gr.GetXaxis().SetLabelOffset(999)
gr.GetYaxis().SetTitle( '(det(I_{ij}) / det(I_{ij}^{full}))^{(1/%s)}'%len(args.variables) if len(args.variables)>1 else 'det(I_{ij}) / det(I_{ij}^{full})' )
gr.GetYaxis().SetLabelSize(0.035)
gr.GetYaxis().SetTitleSize(0.035)
gr.GetYaxis().SetRangeUser(0., 1.1)
gr.Draw( 'AB' )
c1.Update()

# Info
t1 = ROOT.TPaveText(0.55, .85, .95, .92, "blNDC")
t1.SetTextAlign(33)
t1.SetTextSize(0.035)
t1.SetTextFont(1)
t1.SetBorderSize(0)
t1.SetFillColor(ROOT.kWhite)
t1.AddText('Restricted to: ' + ', '.join(args.variables))
t1.AddText('WC: ' + ' '.join(args.parameters))
t1.Draw()

# Selection info
t = ROOT.TLatex()
t.SetTextAngle(90)
t.SetTextSize(0.035)
labels = [ item['plotstring'] for item in data ] 

for i, label in enumerate(labels):
    t.DrawLatex(i+1.25, 0.05, label);

# Directory
plot_directory_ = os.path.join(\
    plot_directory,
    args.plot_directory, 
    sample.name, 
    fisher_directory, 
    args.selection,
    WC_string)

if not os.path.isdir(plot_directory_): os.makedirs(plot_directory_)
c1.Print(os.path.join(plot_directory_, 'fisher_information.png'))

