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
import scipy as spi

# RootTools
from RootTools.core.standard             import *

# TTXPheno
from TTXPheno.Tools.user                 import plot_directory
from TTXPheno.Tools.helpers              import deltaPhi, getCollection, deltaR 
from TTXPheno.Tools.WeightInfo           import WeightInfo

# Import samples
from TTXPheno.samples.benchmarks         import *

# Import helpers
from TTXPheno.Tools.plot_helpers                        import *

# Import process variables
import process_variables

# Import additional
from array                               import array

# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--version',            action='store',      default='v7')
argParser.add_argument('--level',              action='store',      default='gen',  nargs='?', choices=['reco', 'gen', 'genLep'], help='Which level of reconstruction? reco, gen, genLep')
argParser.add_argument('--sample',             action='store',      default='fwlite_ttZ_ll_LO_order2_15weights_ref')
argParser.add_argument('--process',            action='store',      default='ttZ')
argParser.add_argument('--order',              action='store',      default=2)
argParser.add_argument('--small',              action='store_true', help='Run only on a small subset of the data?')
argParser.add_argument('--selection',          action='store',      default='lepSel3-onZ-njet3p-nbjet1p-Zpt0', help="Specify cut.")
argParser.add_argument('--variables',          action='store',      default = [], type=str, nargs='+', help = "argument variables")
argParser.add_argument('--parameters',         action='store',      default = [], type=str, nargs='+', help = "argument parameters")
argParser.add_argument('--luminosity',         action='store',      default=150)

args = argParser.parse_args()

plot_subdirectory = "%s_%s"%(args.level, args.version)
# Import additional functions/classes specified for the level of reconstruction
if args.level == 'reco':     from TTXPheno.Tools.cutInterpreterReco   import cutInterpreter
elif args.level == 'genLep': from TTXPheno.Tools.cutInterpreterGenLep import cutInterpreter
else:                        from TTXPheno.Tools.cutInterpreter       import cutInterpreter

if len(args.parameters) == 0: args.parameters = None

variableColors = [ ROOT.kBlue-7, ROOT.kSpring-5, ROOT.kOrange, ROOT.kPink-1, ROOT.kGray, ROOT.kBlack, ROOT.kYellow ]

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

if len(args.variables) == 0: inputvariables = w.variables
# sort list in same order as w.variables (necessary for EV and EVec calculation)
else: inputvariables = [ item for item in w.variables if item in args.variables ]
# Format input parameters to dict
WC_string = 'SM'
WC = {}
if args.parameters is not None:
    coeffs = args.parameters[::2]
    str_vals = args.parameters[1::2]
    vals   = list( map( float, str_vals ) )
    WC = { coeff:vals[i] for i, coeff in enumerate(coeffs) }
    WC_string = '_'.join( args.parameters )


def get_reweight_function():
    ''' return a weight function
    '''

    def reweight( event, sample ):
        return event.ref_lumiweight1fb * float(args.luminosity) * event_factor

    return reweight


selection_string = cutInterpreter.cutString( args.selection )

# Make sure that weightString contains the same as weightFunction!!!
weightString = 'ref_lumiweight1fb*%s*%s' %( str(args.luminosity) , str(event_factor) )
weightFunction = get_reweight_function()

print weightString
# split selection string step by step
selections = []
selectionElements = args.selection.split('-')
for i, item in enumerate(selectionElements[:-1]):
    selections.append( {'plotstring':' + '.join(replace_selectionstrings(selectionElements[:i+1])) + ' selection', 'selection':'-'.join(selectionElements[:i+1])} )

selections.append( {'plotstring':'full pre-selection (fps)', 'selection':args.selection} )
    
# Additional plot variables from process_variables.py (defined for ttZ, ttW and ttgamma)
plotVariables2D = getattr( process_variables, args.process )['2D']
plotVariables3D = getattr( process_variables, args.process )['3D']

if args.level != 'gen':
    if args.level == 'reco': search_string, replacement_string = ( 'gen', 'reco' )
    # Take care if that is really everything you have to replace!!!
    elif args.level == 'genLep': search_string, replacement_string = ( 'genZ', 'genLepZ' )
    for item in plotVariables2D + plotVariables3D:
        item['var'] = item['var'].replace(search_string, replacement_string)

# Calculate coefficients for binned distribution
# Calculate determinant using the 'variables' submatrix of FI
for var in plotVariables2D:
    #remove initial selection string
    sample.setSelectionString('1')

    var['coeff']       = getCoeffPlotFromDraw( sample, args.order, var['var'], var['binning'], selection_string, weightString=weightString )

    # add bin information to plot labels
    var['plotstring'] = 'fps + ' + var['plotstring'] + ' (%s bins)' %str(var['binning'][0])
    var['color']       = 30

for var in plotVariables3D:
    #remove initial selection string
    sample.setSelectionString('1')

    var['coeff']       = get2DCoeffPlotFromDraw( sample, args.order, var['var'], var['binning'], selection_string, weightString=weightString )
    # add bin information to plot labels
    var['plotstring'] = 'fps + ' + var['plotstring'] + ' (%s:%s bins)' %(str(var['binning'][0]), str(var['binning'][3]))
    var['color']       = 41

# Calculate coefficients unbinned (event loop)
# Calculate determinant using the 'variables' submatrix of FI
for selection in selections:
    selection['coeff'] = getCoeffListFromEvents( sample, selectionString = cutInterpreter.cutString(selection['selection']), weightFunction = weightFunction )
    selection['color']       = 46

# Full Fisher information
full              = { 'plotstring':'full'}
full['coeff']     = getCoeffListFromEvents( sample, selectionString = None, weightFunction = weightFunction )
full['color']     = 15

expo = 1. / len(inputvariables)
data = [full] + selections + plotVariables2D + plotVariables3D
n_data = len(data)

plot_directory_ = os.path.join(\
    plot_directory,
    plot_subdirectory,
    sample.name,
    fisher_directory,
    'eigenvector',
    args.selection,
    WC_string,
    '_'.join(inputvariables) if len([k for k, j in zip(inputvariables, w.variables) if k != j]) > 0 else 'all',
    )
if not os.path.isdir(plot_directory_): os.makedirs(plot_directory_)
if os.path.isfile(os.path.join(plot_directory_,'ev_file.log')): os.remove(os.path.join(plot_directory_,'ev_file.log'))

# Fill dictionaries with Eigenvalues and Eigenvectors
for i,item in enumerate(data):

    print_matrix = w.matrix_to_string( *w.get_total_fisherInformation_matrix( item['coeff'], inputvariables, **WC ) )

    evs, evecs = np.linalg.eigh( w.get_total_fisherInformation_matrix( item['coeff'], inputvariables, **WC )[1] )
    evecs_frac = [ [ abs(entry)/sum(abs(vec)) for entry in vec ] for vec in evecs.T ]
    item['x_pos'] = i+1
    item['evs'] = evs
    item['evecs_frac'] = evecs_frac

    with open(os.path.join(plot_directory_,'ev_file.log'),'a') as f:
        f.write('bin %i\n'%i)
        f.write('matrix\n')
        f.write(print_matrix)
        f.write('\neval\n')
        f.write('(' + '\t'.join([str(item) for item in evs]) + ')')
        f.write('\nevec\n')
        f.write('\n'.join(['(' + '\t'.join([str(el) for el in item]) + ')' for item in evecs.T]))
        f.write('\n\n\n')

# Plots
def drawPlot( log = False ): 
    ''' Plotting function
    '''

    def getTGraph( n, x, y, color=40 ):
        ''' Create a TGraph object
        '''

        gr = ROOT.TGraph( n, x, y )
        gr.SetFillColor( color )
        gr.SetLineWidth( 5 )

        return gr

    # Canvas
    c1 = ROOT.TCanvas()
    c1.SetBottomMargin(0.05)

    mg = ROOT.TGraph( len(data), array( 'd', range( 1, len(data)+1 ) ), array( 'd', [0.]*len(data) ) )
    box = ROOT.TBox( 0, 0, len(data)+1, 0 )
    box.Draw('AL')
    mg.Draw('AB')

    for item in data:
        # create box for each EVal
        for i, y in enumerate(item['evs']):
            xstart = item['x_pos'] - 0.4
            ystart = y * 0.9
            ystop  = y * 1.1

            for j, eVecEntry in enumerate(item['evecs_frac'][i]):
                color = variableColors[j] #eVecEntry[1]
                xstart += (item['evecs_frac'][i][j-1])*0.8 if j>0 else 0
                xstop = xstart + eVecEntry*0.8 #scale to 0.8
                box.SetFillColor(color)
                box.DrawBox(xstart, ystart, xstop, ystop)


    # Layout
    mg.GetXaxis().SetLabelSize(0)
    mg.GetXaxis().SetTickLength(0.)
    mg.GetXaxis().SetLabelOffset(999)
    mg.GetYaxis().SetTitle( 'Eigenvalues I_{ij}' )
    mg.GetYaxis().SetLabelSize(0.035)
    mg.GetYaxis().SetTitleSize(0.035)

    if log: c1.SetLogy()

    # Plot ranges
    ymin = 1e-6
    ymax = 1e6
    mg.GetYaxis().SetRangeUser(ymin, ymax)
    mg.GetXaxis().SetRangeUser(0, len(data)+1)

    c1.Update()

    xstart = 0.95-(len(inputvariables)+1)*0.05-.42

    t1 = ROOT.TPaveText(xstart, .85, .95, .92, "blNDC")
    t1.SetTextAlign(13)
    t1.SetTextSize(0.035)
    t1.SetTextFont(1)
    t1.SetBorderSize(0)
    t1.SetFillColor(ROOT.kWhite)
    t1.AddText('Eigenvector composition: ')
    t1.AddText('WC: ' + ' '.join(args.parameters) if args.parameters is not None else 'Standard Model' )
    t1.Draw()

    xstart += 0.42
    text = []
    # Info
    for i, var in enumerate(inputvariables):
        if i != 0: xstart += (len(inputvariables[i-1])+2)*0.016
        ti = ROOT.TPaveText(xstart, .85, 1, .92, "blNDC")
        ti.SetTextAlign(13)
        ti.SetTextSize(0.035)
        ti.SetTextFont(1)
        ti.SetBorderSize(0)
        ti.SetFillColor(ROOT.kWhite)
        ti.SetFillStyle(3050)
        ti.AddText(var)
        ti.GetListOfLines().Last().SetTextColor(variableColors[i])
        ti.AddText('')
        text.append(ti)
        del ti

    for item in text:
        item.Draw()

    # Selection info
    t = ROOT.TLatex()
    t.SetTextAngle(90)
    t.SetTextSize(0.035)

    # Labeling
    for i, item in enumerate(data):
        t.DrawLatex( i+1.25, ymin*1.4 if log else 0.05, "#color[%i]{%s}"%(item['color'],item['plotstring']) )

    # Directory
    plot_directory_ = os.path.join(\
        plot_directory,
        plot_subdirectory,
        sample.name,
        fisher_directory,
        'eigenvector',
        args.selection,
        WC_string,
        '_'.join(inputvariables) if len([i for i, j in zip(inputvariables, w.variables) if i != j]) > 0 else 'all',
        'log' if log else 'lin')
    
    if not os.path.isdir(plot_directory_): os.makedirs(plot_directory_)
    c1.Print(os.path.join(plot_directory_, 'eigenvector.png'))
    
    del c1


# Plot lin and log plots
for log in [ True, False ]:
    drawPlot( log )

