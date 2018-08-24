''' Fit EFT parameters to observation.
    Currently, this class is a placeholder to interface to Wolfgangs code.
'''

# Standard imports 
import sys, copy
import ROOT
import imp
import pickle
import gc

import numpy as np

from multiprocessing import Pool

from math import sqrt
# turn off graphics
ROOT.gROOT.SetBatch( True )

# RootTools
from RootTools.core.standard import *

# Logger
import TTXPheno.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   'DEBUG', logFile = None)
logger_rt = logger_rt.get_logger('INFO', logFile = None)

# TTXPheno
from TTXPheno.samples.benchmarks import * 
from TTXPheno.Tools.user import plot_directory

# get the reweighting function
from TTXPheno.Tools.WeightInfo import WeightInfo

# Load the analysis regions
#from TTXPheno.Analysis.regions import ttZRegions as regions
from TTXPheno.Analysis.regions import ttgammaRegions as regions

# Arguments
import argparse

argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--version',         action='store',     default='test', help='Appendix to plot directory')
argParser.add_argument('--process',         action='store',     default='ttgamma_1l', nargs='?', choices=['ttZ_3l', 'ttZ_4l', 'ttgamma_1l', 'ttgamma_2l'], help="which process to calculate?")
argParser.add_argument('--sample',          action='store',     default='fwlite_ttgamma_LO_order2_15weights_ref', help='Sample name specified in sample/python/benchmarks.py, e.g. fwlite_ttZ_ll_LO_order2_15weights_ref')
argParser.add_argument('--order',           action='store',     default=2, help='Polynomial order of weight string (e.g. 2)')
argParser.add_argument('--selection',       action='store',     default='lepSel1-gammapt40-njet3p-nbjet1p-relIso0to0.4-met40', help="Specify cut.")
argParser.add_argument('--small',           action='store_true', help='Run only on a small subset of the data?')
argParser.add_argument('--level',           action='store',     default='gen', nargs='?', choices=['reco', 'gen'], help='Which level of reconstruction? reco, gen')
argParser.add_argument('--variables' ,      action='store',     default = ['ctZ', 'ctZI'], type=str, nargs='+', help = "argument plotting variables")
argParser.add_argument('--binning',         action='store',     default = [21, -10, 10, 21, -10, 10], type=int, nargs='+', help = "argument parameters")
argParser.add_argument('--luminosity',      action='store',     default=150, help='Luminosity for weighting the plots')

args = argParser.parse_args()

if len(args.binning) != 6:
    raise ValueError('Binning must be 6 values in the form BINSX STARTX STOPX BINSY STARTY STOPY! Input: %s' %' '.join(args.binning))

#if (binningX[2] - binningX[1]) % binningX[0] != 0:
#    raise ValueError('Binning: Difference of upper and lower bound must be a multiple of the number of bins: (%s - %s) / %s != integer'%(str(binningX[2]), str(binningX[1]), str(binningX[0])) )

#if (binningY[2] - binningY[1]) % binningY[0] != 0:
#    raise ValueError('Binning: Difference of upper and lower bound must be a multiple of the number of bins: (%s - %s) / %s != integer'%(str(binningY[2]), str(binningY[1]), str(binningY[0])) )

if len(args.variables) != 2:
    raise ValueError('Give TWO input variables in --variables! Input: %s' %' '.join(args.variables))

# Import additional functions/classes specified for the level of reconstruction
if args.level == 'reco':
    from TTXPheno.Tools.cutInterpreterReco import cutInterpreter
else:
    from TTXPheno.Tools.cutInterpreterGen import cutInterpreter


if args.process == 'ttgamma_1l': ttSampleName = 'fwlite_tt_nonhad_LO_order2_15weights'
else: ttSampleName = 'fwlite_tt_dilep_LO_order2_15weights'

# Import samples
sample_file = "$CMSSW_BASE/python/TTXPheno/samples/benchmarks.py"
loadedSamples = imp.load_source( "samples", os.path.expandvars( sample_file ) )

ttXSample = getattr( loadedSamples, args.sample )
WZSample        = getattr( loadedSamples, 'fwlite_WZ_lep_LO_order2_15weights' )
ttSample        = getattr( loadedSamples, ttSampleName )
tWSample        = getattr( loadedSamples, 'fwlite_tW_LO_order2_15weights' )
tWZSample       = getattr( loadedSamples, 'fwlite_tWZ_LO_order2_15weights' )
tZqSample       = getattr( loadedSamples, 'fwlite_tZq_LO_order2_15weights' )
ZgammaSample    = getattr( loadedSamples, 'fwlite_Zgamma_LO_order2_15weights' )
ttgammaSample   = getattr( loadedSamples, 'fwlite_ttgamma_bg_LO_order2_15weights' )


if args.process.split('_')[0] == 'ttgamma':
    ttgammaIsrSample  = copy.deepcopy( ttXSample ) #select ttgamma events with isolated gamma from ISR (cat a2)
    ttgammaIsrSample.name = 'fwlite_ttgamma_ISR_LO_order2_15weights_ref'

elif args.process.split('_')[0] == 'ttZ':
    ttZIsrSample  = copy.deepcopy( ttXSample ) #Z not from t/gluons
    ttZIsrSample.name = 'fwlite_ttZ_ISR_LO_order2_15weights_ref'

if args.process == 'ttZ_3l': bg = [ ttZIsrSample, WZSample, tWZSample, tZqSample, ttgammaSample ]
elif args.process == 'ttZ_4l': bg = [ ttZIsrSample, WZSample, tWZSample, tZqSample, ttgammaSample ]
elif args.process == 'ttgamma_1l': bg = [ ttSample, ttgammaIsrSample, tWSample, tWZSample, tZqSample, ZgammaSample ]
elif args.process == 'ttgamma_2l': bg = [ ttSample, ttgammaIsrSample, tWSample, tWZSample, tZqSample, ZgammaSample ]

if args.small:
    for s in [ttXSample] + bg:
        s.reduceFiles( to = 20 )

def checkReferencePoint( sample ):
    ''' check if sample is simulated with a reference point
    '''
    return pickle.load(file(sample.reweight_pkl))['ref_point'] != {}

# set selection string
selectionString = cutInterpreter.cutString(args.selection)

# configure samples
for s in [ttXSample] + bg:

    s.event_factor = s.nEvents / float( s.chain.GetEntries() )
    s.weightInfo = WeightInfo( s.reweight_pkl )
    s.weightInfo.set_order( args.order )
    s.setSelectionString( selectionString )

    if checkReferencePoint( s ):
#        s.read_variables = ["ref_lumiweight1fb/F", VectorTreeVariable.fromString('p[C/F]', nMax=2000) ]
        s.setWeightString( 'ref_lumiweight1fb*(%s)*(%s)'%( str(args.luminosity), str(s.event_factor) ) )
#        s.setWeightString( '(%s)*1000/(%s)/p_C[0]*(%s)*(%s)'%( str(s.xsec), str(s.nEvents), str(args.luminosity), str(s.event_factor) ) )
    else:
#        s.read_variables = ["lumiweight1fb/F"]
        s.setWeightString( 'lumiweight1fb*(%s)*(%s)'%( str(args.luminosity), str(s.event_factor) ) )
#        s.setWeightString( '(%s)*1000/(%s)*(%s)*(%s)'%( str(s.xsec), str(s.nEvents), str(args.luminosity), str(s.event_factor) ) )

# overlap removal and signal selection
if args.process.split('_')[0] == 'ttgamma':
    ttXSample.addSelectionString(        "(signalPhoton==1)"               ) #cat a1 <- the one and only signal
    ttgammaIsrSample.addSelectionString( "(isrPhoton==1)"                  ) #cat a2
    ttSample.addSelectionString(         "(signalPhoton!=1&&isrPhoton!=1)" ) #cat b,c,d

elif args.process.split('_')[0] == 'ttZ':
    ttXSample.addSelectionString(    "(signalZ==1)"               ) 
    ttZIsrSample.addSelectionString( "(signalZ==0)"                    ) 


for var in args.variables:
    if var not in ttXSample.weightInfo.variables:
        raise ValueError('Input variable not in gridpack: %s' %var)

observation                  = {}
signal_jec_uncertainty       = {}
signal_fakerate_uncertainty  = {}

ttX_SM_rate                  = {}
ttX_SM_jec_uncertainty       = {}
ttX_SM_fakerate_uncertainty  = {}
ttX_coeffList                = {}

background_rate                 = {}
background_jec_uncertainty      = {}
background_fakerate_uncertainty = {}

for i_region, region in enumerate(regions):
    # compute signal yield for this region (this is the final code)

    logger.info( "At region %s", region )

    # ttX SM
    ttX_coeffList[region] = ttXSample.weightInfo.getCoeffListFromDraw( ttXSample, selectionString = region.cutString() )
    ttX_SM_rate[region]   = ttXSample.weightInfo.get_weight_yield( ttX_coeffList[region] )

    ttX_SM_jec_uncertainty      [region] = 1.05 
    ttX_SM_fakerate_uncertainty [region] = 1.0  # signal has no FR uncertainty

    # signal uncertainties
    signal_jec_uncertainty      [region] = 1.05 
    signal_fakerate_uncertainty [region] = 1.0  # signal has no FR uncertainty

    background_rate[region]                 = {}
    background_fakerate_uncertainty[region] = {}
    background_jec_uncertainty[region]      = {}

    for i_background, background in enumerate(bg):
        # compute bg yield for this region (this is the final code)

        background_rate                 [region][background.name] = background.getYieldFromDraw( selectionString=region.cutString() )['val']
        background_fakerate_uncertainty [region][background.name] = 1   + 0.03*i_region*(i_background+1) #change that
        background_jec_uncertainty      [region][background.name] = 1.2 - 0.02*i_region*(i_background+1) #change that

    # Our expected observation :-)
    observation[region] = int( sum( background_rate[region].values() ) + ttX_SM_rate[region] )

# Write temporary card file
from TTXPheno.Tools.cardFileWriter import cardFileWriter
c = cardFileWriter.cardFileWriter()

def calculation( variables ):
#def calculation( var1, var2 ):

        var1, var2 = variables
        kwargs = { args.variables[0]:var1, args.variables[1]:var2 }

        # uncertainties
        c.reset()
        c.addUncertainty('lumi', 'lnN')
        c.addUncertainty('JEC',  'lnN')
        c.addUncertainty('fake', 'lnN')

        #processes = [ "signal", "ttX_SM_yield" ] + [s.name for s in bg]
        signal_rate                  = {}
        for i_region, region in enumerate(regions):
 
            signal_rate[region] = ttXSample.weightInfo.get_weight_yield( ttX_coeffList[region], **kwargs) - ttX_SM_rate[region] 

            bin_name = "Region_%i" % i_region
            nice_name = region.__str__()
            c.addBin(bin_name,['ttX_SM'] + ['_'.join(s.name.split('_')[1:3]) for s in bg], nice_name)
            c.specifyObservation( bin_name, observation[region] )

            c.specifyFlatUncertainty( 'lumi', 1.05 )

            c.specifyExpectation( bin_name, 'signal', signal_rate[region] )
            c.specifyUncertainty( 'JEC', bin_name, 'signal', signal_jec_uncertainty[region])
            c.specifyUncertainty( 'fake',bin_name, 'signal', signal_fakerate_uncertainty[region])

            c.specifyExpectation( bin_name, 'ttX_SM', ttX_SM_rate[region] )
            c.specifyUncertainty( 'JEC', bin_name, 'ttX_SM', ttX_SM_jec_uncertainty[region])
            c.specifyUncertainty( 'fake',bin_name, 'ttX_SM', ttX_SM_fakerate_uncertainty[region])

            for background in bg:
                c.specifyExpectation( bin_name, '_'.join(background.name.split('_')[1:3]), background_rate[region][background.name] )
                c.specifyUncertainty( 'JEC', bin_name, '_'.join(background.name.split('_')[1:3]), background_jec_uncertainty[region][background.name])
                c.specifyUncertainty( 'fake',bin_name, '_'.join(background.name.split('_')[1:3]), background_fakerate_uncertainty[region][background.name])
                
        nameList = ttXSample.name.split('_')[1:3] + args.variables + args.binning + [ args.level, args.version, args.order, args.luminosity, args.selection, 'small' if args.small else 'full', var1, var2,  ]
        cardname = '%s_limit_card'%'_'.join( map( str, nameList ) )
        c.writeToFile( './tmp/%s.txt'%cardname )

        # try to adjust rmax with some margin
        exp_tot_sigmas = 0
        max_rmax = float('inf')
        for region in regions:

            tot_background = sum( [ background_rate[region][background.name] for background in bg ] )
            exp_tot_sigmas += abs(signal_rate[region]) / sqrt( tot_background ) if tot_background > 0 else 0

            # avoid total neg. yield
            if signal_rate[region] < 0:
                max_r = -tot_background / signal_rate[region]
                if max_r < max_rmax:
                    max_rmax = max_r
            
        rmax_est = 400. / exp_tot_sigmas if exp_tot_sigmas > 0 else 1.
        if max_rmax < rmax_est:
            rmax_est = 0.9*max_rmax # safety margin such that at least +10% total yield survives in the smallest SR

        profiledLoglikelihoodFit = ProfiledLoglikelihoodFit( './tmp/%s.txt'%cardname )
        profiledLoglikelihoodFit.make_workspace(rmin=0, rmax=rmax_est)
        #expected_limit = profiledLoglikelihoodFit.calculate_limit( calculator = "frequentist" )
        expected_limit = profiledLoglikelihoodFit.calculate_limit( calculator = "asymptotic" )
        logger.info( "Expected Limit: %f", expected_limit[0] )
        profiledLoglikelihoodFit.cleanup()
        del profiledLoglikelihoodFit
        ROOT.gDirectory.Clear()

        return var1, var2, expected_limit[0]



# Limit plot
from TTXPheno.Analysis.ProfiledLoglikelihoodFit import ProfiledLoglikelihoodFit

binningX = args.binning[:3]
binningY = args.binning[3:]

xRange = np.linspace( binningX[1], binningX[2], binningX[0], endpoint=False)
yRange = np.linspace( binningY[1], binningY[2], binningY[0], endpoint=False)

xRange = [ el+ 0.5*(xRange[1]-xRange[0]) for el in xRange ]
yRange = [ el+ 0.5*(yRange[1]-yRange[0]) for el in yRange ]

results = []
for varX in xRange:
    # do not run all calc in one pool, memory leak!!!
    pool = Pool( processes = 8 )
    results += pool.map( calculation, [ (varX, varY) for varY in yRange ] )
    del pool

limit = ROOT.TH2F( 'limit', 'limit', *(binningX + binningY) )
limit.GetZaxis().SetTitle("-2 #Delta ln(L)")
#limit.GetZaxis().SetTitleOffset(0.8)
limit.GetZaxis().SetLabelSize(0.045)
limit.GetZaxis().SetTitleSize(0.045)



#scale to SM
results = [ (el[0], el[1], np.log(el[2])) for el in results ]
results.sort( key = lambda res: -res[2] )
print results
maxResult = results[0][2]

#results.sort( key = lambda res: (res[0],res[1],-res[2]) )
#use the point closesed to the SM as reference
#maxResult = results[1][2] #do not use SM point (pole at 0,0)

for x,y,result in results:
    limit.SetBinContent( limit.FindBin(x, y), -2*(result - maxResult) )

limitPlot = Plot2D.fromHisto( '_'.join(args.variables), texX = args.variables[0], texY = args.variables[1], histos = [[limit]] )
limitPlot.drawOption = "colz"
ROOT.gStyle.SetPaintTextFormat("2.2f")

def drawObjects( hasData = False ):
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextSize(0.04)
    tex.SetTextAlign(11) # align right
    lines = [
      (0.15, 0.95, "Simulation (gen)"),
      (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV)' %float(args.luminosity) )
    ]
    return [tex.DrawLatex(*l) for l in lines]

plot_directory_ = os.path.join(\
    plot_directory,
    '%s_%s'%(args.level, args.version),
    ttXSample.name,
    'backgrounds',
    'limit_small' if args.small else 'limit',
    args.selection)

plotting.draw2D( limitPlot,
                 logZ=False,
                 plot_directory = plot_directory_,
                 drawObjects = drawObjects(),
                 extensions = ["png"])

