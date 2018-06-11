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
argParser.add_argument('--selection',          action='store',      default='lepSel3-onZ-njet3p-nbjet1p-Zpt0', help="Specify cut.")
argParser.add_argument('--small',              action='store_true', help='Run only on a small subset of the data?')

argParser.add_argument('--parameters',         action='store',      default = [], type=str, nargs='+', help = "argument parameters")
#argParser.add_argument('--parameters',         action='store',      default = ['ctW', '3', 'ctWI', '3', 'ctZ', '3', 'ctZI', '3'], type=str, nargs='+', help = "argument parameters")
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
weightString = 'ref_lumiweight1fb*%s*%s' %( str(args.luminosity), str(event_factor) )


coeffs = []
str_vals = []
vals = []
WC = {}
WC_string = 'SM'
if len(args.parameters) != 0:
    coeffs = args.parameters[::2]
    str_vals = args.parameters[1::2]
    vals   = list( map( float, str_vals ) )
    WC = { coeff:vals[i] for i, coeff in enumerate(coeffs) }
    WC_string = '_'.join(args.parameters )
    weightString += '*' + w.get_weight_string( **WC )
    print(weightString)

variables = [
             'cpQM',
             'cpt',
             'ctZ',
             'ctZI'
]


plot_variables = [
                  { 'plotstring':'p_{T}(Z) + sel',       'name':'Z_pt',           'binning':[ 20, 0, 400 ] },
                  { 'plotstring':'cos(theta^*) + sel',  'name':'Z_cosThetaStar', 'binning':[ 20, -1.2, 1.2 ] },
                  { 'plotstring':'m_{ll} + sel',         'name':'Z_mass',         'binning':[ 20, 70, 110 ] },
                  { 'plotstring':'phi(Z) + sel',        'name':'Z_phi',          'binning':[ 20, -np.pi, np.pi ] },
                  { 'plotstring':'eta(Z) + sel',        'name':'Z_eta',          'binning':[ 20, -3, 3 ] },
]


selectionStrings = [
                  { 'plotstring':'sel',         'selection':'lepSel3-onZ-njet3p-nbjet1p-Zpt0' },
                  { 'plotstring':'sel 0-100',   'selection':'lepSel3-onZ-njet3p-nbjet1p-Zpt0to100' },
                  { 'plotstring':'sel 100-200', 'selection':'lepSel3-onZ-njet3p-nbjet1p-Zpt100to200' },
                  { 'plotstring':'sel 200-300', 'selection':'lepSel3-onZ-njet3p-nbjet1p-Zpt200to300' },
                  { 'plotstring':'sel 300-400', 'selection':'lepSel3-onZ-njet3p-nbjet1p-Zpt300to400' },
                  { 'plotstring':'sel 400',     'selection':'lepSel3-onZ-njet3p-nbjet1p-Zpt400' },
                  { 'plotstring':'3 leps',      'selection':'lepSel3' },
                  { 'plotstring':'onZ',         'selection':'onZ' },
                  { 'plotstring':'3 jets',      'selection':'njet3p' },
                  { 'plotstring':'1 bjet',      'selection':'nbjet1p' },
]


coeff_var =  [ getCoeffPlotFromDraw( sample, args.order, var['name'], var['binning'], selection_string, weightString=weightString ) for var in plot_variables ]
coeff_sels = [ getCoeffListFromEvents( sample, cutInterpreter.cutString(selection['selection']), args.luminosity*event_factor ) for selection in selectionStrings ]
coeff_tot =  getCoeffListFromEvents( sample, None, args.luminosity*event_factor )
#coeff_sel = getCoeffListFromEvents( sample, selection_string, args.luminosity*event_factor )

detI_var =  [ np.linalg.det( w.get_total_fisherInformation_matrix( coeff, variables )[1] ) for coeff in coeff_var ]
detI_sels = [ np.linalg.det( w.get_total_fisherInformation_matrix( c_sel, variables )[1] ) for c_sel in coeff_sels ]
detI_tot =  np.linalg.det( w.get_total_fisherInformation_matrix( coeff_tot, variables )[1] )
#detI_sel = np.linalg.det( w.get_total_fisherInformation_matrix( coeff_sel, variables )[1] )

det = np.array([detI_tot] + detI_sels +  detI_var) / detI_tot
n = len(det)

print (det)

from array import array
x = array('d')
y = array('d')

for i in range(n):
    x.append(i+1)
    y.append( (det[i])**(1./len(variables)) )

c1 = ROOT.TCanvas()
gr = ROOT.TGraph( n, x, y )
gr.SetFillColor( 40 )
gr.SetLineWidth( 4 )
gr.GetXaxis().SetLabelSize(0)
gr.GetXaxis().SetLabelOffset(999)
#gr.SetMarkerColor( 4 )
#gr.SetMarkerStyle( 21 )
#gr.SetTitle( 'a simple graph' )
#gr.GetXaxis().SetTitle( 'X title' )
gr.GetYaxis().SetTitle( '(det(I_{ij}) / det(I_{ij}^{full}))^{(1/%s)}'%len(variables) if len(variables)>1 else 'det(I_{ij}) / det(I_{ij}^{full})' )
gr.GetYaxis().SetLabelSize(0.035)
gr.GetYaxis().SetTitleSize(0.035)
gr.Draw( 'AB' )
c1.Update()


t1 = ROOT.TText()
t1.SetTextSize(0.035)
label = 'Restricted to: ' + ', '.join(variables)
t1.DrawText( len(selectionStrings + plot_variables) - len(variables), 1., label);

t2 = ROOT.TText()
t2.SetTextSize(0.035)
label = ' '.join(args.parameters)
t2.DrawText( len(selectionStrings + plot_variables) - .3*len(args.parameters), .95, label);


t = ROOT.TText()
t.SetTextAngle(90)
#t.SetTextAlign(32)
t.SetTextSize(0.035)
#t.SetTextFont(72)
labels = ['full'] + [item['plotstring'] for item in selectionStrings + plot_variables]

for i, label in enumerate(labels):
    t.DrawText(i+1.25, 0.2, label);


plot_directory_ = os.path.join(\
    plot_directory,
    args.plot_directory, 
    sample.name, 
    'fisher_information', 
    args.selection,
    WC_string)

if not os.path.isdir(plot_directory_): os.makedirs(plot_directory_)
c1.Print(os.path.join(plot_directory_, 'fisher_information.png'))

