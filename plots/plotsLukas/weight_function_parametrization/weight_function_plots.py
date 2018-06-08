''' simple draw script
'''

# Standrard imports
import ROOT
import os, imp

# numpy
import numpy as np

# RootTools
from RootTools.core.standard import *

# Logger
import TTXPheno.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger = logger.get_logger( 'INFO', logFile = None )

# TTXPheno
from TTXPheno.Tools.user           import plot_directory
from TTXPheno.Tools.WeightInfo     import WeightInfo
from TTXPheno.Tools.cutInterpreter import cutInterpreter

# Sample
from TTXPheno.samples.benchmarks   import *
from plot_helpers                  import *

import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--plot_directory',     action='store',      default='weight_function')
argParser.add_argument('--sample',             action='store',      default='fwlite_ttZ_ll_LO_order2_15weights_ref')
argParser.add_argument('--order',              action='store',      default=2)
argParser.add_argument('--selection',          action='store',      default='lepSel3-onZ-njet3p-nbjet1p-Zpt0', help="Specify cut.")
argParser.add_argument('--small',              action='store_true', help='Run only on a small subset of the data?')
argParser.add_argument('--differential',       action='store_true', help='Create differential plots')
argParser.add_argument('--normalized',         action='store_true', help='Normalize plots to SM')
argParser.add_argument('--parameters',         action='store',      default = ['ctW', 'ctWI', 'ctZ', 'ctZI'], type=str, nargs='+', help = "argument parameters")

args = argParser.parse_args()

if args.differential and args.normalized: exit()
textsize = 0.04

# Import samples
sample_file = "$CMSSW_BASE/python/TTXPheno/samples/benchmarks.py"
samples     = imp.load_source( "samples", os.path.expandvars( sample_file ) )
sample      = getattr( samples, args.sample )

if args.small: sample.reduceFiles( to = 1 )

# Make subdirectory
subDirectory = []
if args.small:     subDirectory.append("small")
if args.differential: subDirectory.append("differential")
if args.normalized: subDirectory.append("normalized")
subDirectory = '_'.join( subDirectory )

# Plot directory
plot_directory_ = os.path.join(\
    plot_directory,
    args.plot_directory,
    sample.name,
    subDirectory)

if not os.path.isdir(plot_directory_): os.makedirs(plot_directory_)

# get TChain
#c = sample.chain 

# Polynomial parametrization
w = WeightInfo( sample.reweight_pkl )
w.set_order( int(args.order) )

colors = [ ROOT.kMagenta+1, ROOT.kOrange, ROOT.kBlue, ROOT.kCyan+1, ROOT.kGreen+1, ROOT.kRed, ROOT.kViolet, ROOT.kYellow+2 ]

# selection
selection_string = cutInterpreter.cutString(args.selection)

# sample weight
#weight_string = '150*lumiweight1fb'
weight_string = None

# coefficient list (all weights)
#coeffList =  getCoeffPlotFromDraw( sample, args.order, 'Z_pt', [ 20, 0, 500 ], selection_string, weightString=weight_string)
coeffList = [getCoeffListFromDraw( sample, args.order, selection_string, weightString = weight_string )]

def plot_1D_yield( var, coeffList, range = [-10, 10, .5], differential = False, normalization = False ):
    ''' Create 1D plot of (differential) yield as a function of WC (val) for a given list of weights (coeffList)
    '''

    var_val = np.arange( range[0], range[1], range[2] )
    dict_list = [ {var:val} for val in var_val]
    if differential:
        y = np.array( [ w.get_diff_weight_yield( var, coeffList, **item ) for item in dict_list] )
    else:
        y = np.array( [ w.get_weight_yield( coeffList, **item ) for item in dict_list] )
        if normalization: y /= w.get_weight_yield( coeffList, **{var:0} )

    graph = ROOT.TGraph(len(var_val), var_val, y)

    graph.SetLineWidth(2)
    graph.SetLineColor(ROOT.kRed)
    graph.GetXaxis().SetRangeUser(range[0],range[1])

    y_label = ''
    if differential: y_label += 'Diff. '
    y_label += 'Total Yield'
    if normalization and not differential: y_label += ' [SM Yield]'

    graph.GetYaxis().SetTitle(y_label)

    graph.GetXaxis().SetTitle(var)
    graph.SetTitle('')

    graph.GetXaxis().SetTitleSize(textsize)
    graph.GetYaxis().SetTitleSize(textsize)
    graph.GetXaxis().SetLabelSize(textsize)
    graph.GetYaxis().SetLabelSize(textsize)

    # Base Line
    if normalization:
        line = ROOT.TLine(range[0],1,range[1],1)
    if differential:
        line = ROOT.TLine(range[0],0,range[1],0)
    if normalization or differential:
        line.SetLineColor(ROOT.kBlack)
        line.SetLineWidth(2)
        line.SetLineColor(ROOT.kBlue)

    # Legend
    legend = ROOT.TLegend(0.63, 0.78, .9, .9)
    legend.SetTextSize(textsize)
    legend.SetTextFont(1)
    if normalization: legend.AddEntry(line, "SM Yield", "l")
    legend.AddEntry(graph, "Total Yield", "l")

    # Plotting
    c1 = ROOT.TCanvas()
    graph.Draw()
    if normalization or differential: line.Draw('SAME')
    legend.Draw()
    logger.info('Plot created:' + os.path.join(plot_directory_, '%s.png') %(var))
    c1.Print(os.path.join(plot_directory_, '%s.png') %(var))


for var in args.parameters:
    plot_1D_yield(var, coeffList, range = [-10, 10, .1], differential = args.differential, normalization = args.normalized)
