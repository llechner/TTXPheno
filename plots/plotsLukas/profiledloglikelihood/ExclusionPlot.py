''' Fit EFT parameters to observation.
    Currently, this class is a placeholder to interface to Wolfgangs code.
'''

# Standard imports 
import sys
import ROOT
import imp
import pickle

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
from TTXPheno.Analysis.regions import regions

# Arguments
import argparse

argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--version',         action='store',     default='v26', help='Appendix to plot directory')
argParser.add_argument('--sample',          action='store',     default='fwlite_ttZ_ll_LO_order2_15weights_ref', help='Sample name specified in sample/python/benchmarks.py, e.g. fwlite_ttZ_ll_LO_order2_15weights_ref')
argParser.add_argument('--order',           action='store',     default=2, help='Polynomial order of weight string (e.g. 2)')
argParser.add_argument('--selection',       action='store',     default='lepSel3-onZ-njet3p-nbjet1p-Zpt0', help="Specify cut.")
argParser.add_argument('--small',           action='store_true', help='Run only on a small subset of the data?')
argParser.add_argument('--parameters',      action='store',     default = [], type=str, nargs='+', help = "argument parameters")
argParser.add_argument('--level',           action='store',     default='gen', nargs='?', choices=['reco', 'gen', 'genLep'], help='Which level of reconstruction? reco, gen, genLep')
argParser.add_argument('--variables' ,      action='store',     default = ['cpQM', 'cpt'], type=str, nargs='+', help = "argument plotting variables")
argParser.add_argument('--binning',         action='store',     default = [18, -8, 28, 18, -23, 13], type=int, nargs='+', help = "argument parameters")
argParser.add_argument('--luminosity',      action='store',     default=150, help='Luminosity for weighting the plots')

args = argParser.parse_args()

if len(args.binning) != 6:
    raise ValueError('Binning must be 6 values in the form BINSX STARTX STOPX BINSY STARTY STOPY! Input: %s' %' '.join(args.binning))

if len(args.variables) != 2:
    raise ValueError('Give TWO input variables in --variables! Input: %s' %' '.join(args.variables))

# Import additional functions/classes specified for the level of reconstruction
if args.level == 'reco':
    from TTXPheno.Tools.cutInterpreterReco import cutInterpreter
elif args.level == 'genLep':
    from TTXPheno.Tools.cutInterpreterGenLep import cutInterpreter
else:
    from TTXPheno.Tools.cutInterpreterGen import cutInterpreter

# Import samples
sample_file = "$CMSSW_BASE/python/TTXPheno/samples/benchmarks.py"
loadedSamples = imp.load_source( "samples", os.path.expandvars( sample_file ) )
ttXSample = getattr( loadedSamples, args.sample )
WZSample = getattr( loadedSamples, 'fwlite_WZ_lep_LO_order2_15weights' )
ttSample = getattr( loadedSamples, 'fwlite_tt_lep_LO_order2_15weights' )
ttSemiLepSample = getattr( loadedSamples, 'fwlite_tt_semilep_LO_order2_15weights' )
tWSample = getattr( loadedSamples, 'fwlite_tW_LO_order2_15weights' )
tWZSample = getattr( loadedSamples, 'fwlite_tWZ_LO_order2_15weights' )
tZqSample = getattr( loadedSamples, 'fwlite_tZq_LO_order2_15weights' )
ZgammaSample = getattr( loadedSamples, 'fwlite_Zgamma_LO_order2_15weights' )
ttgammaSample = getattr( loadedSamples, 'fwlite_ttgamma_bg_LO_order2_15weights' )

# Polynomial parametrization
# ATTENTION IF U USE MORE THAN ONE SIGNAL SAMPLE!!!
w = WeightInfo(ttXSample.reweight_pkl)
w.set_order(int(args.order))

for var in args.variables:
    if var not in w.variables:
        raise ValueError('Input variable not in gridpack: %s' %var)

# set selection string
selectionString = cutInterpreter.cutString(args.selection)

signal = ttXSample
bg = [ WZSample, ttSample ]

if args.small:
    for s in [signal] + bg:
        s.reduceFiles( to = 5 )

def checkReferencePoint( sample ):
    ''' check if sample is simulated with a reference point
    '''
    return pickle.load(file(sample.reweight_pkl))['ref_point'] != {}

# configure samples
for s in [signal] + bg:
    # calculate compensation factor (if sample is reduce due to small or some file weren't processed)
    s.event_factor = s.nEvents / float( s.chain.GetEntries() )
    s.setSelectionString( selectionString )
    if checkReferencePoint( s ):
        s.read_variables = ["ref_lumiweight1fb/F"]
        s.read_variables.append( VectorTreeVariable.fromString('p[C/F]', nMax=2000) )
        s.setWeightString( 'ref_lumiweight1fb*(%s)*(%s)*(%s)'%(str(args.luminosity), str(s.event_factor), w.get_weight_string(**{})) )
    else:
        s.read_variables = ["lumiweight1fb/F"]
        s.setWeightString( 'lumiweight1fb*(%s)*(%s)'%(str(args.luminosity), str(s.event_factor)) )

# REMOVE THAT WHEN BG IMPLEMENTED!!!
class FakeBackground:
    ''' Just a fake background number producer. Ignores everything, returns a number.
    '''
    def __init__( self, name, relative_fraction):
        self.name              = name
        self.relative_fraction = relative_fraction
    def setSelectionString( self, *args, **kwargs):
        pass 

# Backgrounds
bg = [FakeBackground("bkg1",.1), FakeBackground("bkg2", .2)]

# parameter point
if len(args.parameters) < 2: args.parameters = None

params = {}
if args.parameters is not None:
    coeffs = args.parameters[::2]
    str_vals = args.parameters[1::2]
    vals = list( map( float, str_vals ) )

    for (coeff, val, str_val) in zip(coeffs, vals, str_vals):
        if coeff not in w.variables:
            raise ValueError('Parameter not in Gridpack: %s'%coeff)
        params[ coeff ] = val


observation                  = {}
signal_jec_uncertainty       = {}
signal_fakerate_uncertainty  = {}
ttX_SM_rate                  = {}
ttX_SM_jec_uncertainty       = {}
ttX_SM_fakerate_uncertainty  = {}
background_rate                 = {}
background_jec_uncertainty      = {}
background_fakerate_uncertainty = {}

ttX_coeffList                = {}
for i_region, region in enumerate(regions):

    logger.info( "At region %s", region )

    # compute signal yield for this region (this is the final code)
    ttX_coeffList[region] = w.getCoeffListFromDraw( signal, region.cutString(), weightString=signal.weightString )
    # ttX SM
    ttX_SM_rate[region] = w.get_weight_yield( ttX_coeffList[region] )
    ttX_SM_jec_uncertainty      [region] = 1.05 
    ttX_SM_fakerate_uncertainty [region] = 1.0  # signal has no FR uncertainty

    # signal uncertainties
    signal_jec_uncertainty      [region] = 1.05 
    signal_fakerate_uncertainty [region] = 1.0  # signal has no FR uncertainty

    background_rate[region]                = {}
    background_fakerate_uncertainty[region] = {}
    background_jec_uncertainty[region]      = {}
    for i_background, background in enumerate(bg):
        # compute some  background yield
        background_rate[region][background.name] = background.relative_fraction*ttX_SM_rate[region]

        # Guess the uncertainty
        background_fakerate_uncertainty[region][background.name] =  1   + 0.03*i_region*(i_background+1) 
        background_jec_uncertainty[region][background.name]      =  1.2 - 0.02*i_region*(i_background+1) 

    # Our expected observation :-)
    observation[region] = int( sum( background_rate[region].values() ) + ttX_SM_rate[region] )


# Write temporary card file
from TTXPheno.Tools.cardFileWriter import cardFileWriter
c = cardFileWriter.cardFileWriter()

# Limit plot
from TTXPheno.Analysis.ProfiledLoglikelihoodFit import ProfiledLoglikelihoodFit

#binningX = binningY = [12, -12, 12]
binningX = args.binning[:3] #[24, -24-6, 24-6]
binningY = args.binning[3:] #[24, -24+16, 24+16]

limit = ROOT.TH2F( 'limit', 'limit', *(binningX + binningY) )

for var1 in range( binningX[1], binningX[2], ( binningX[2] - binningX[1]) / binningX[0] ):
    for var2 in range( binningY[1], binningY[2], ( binningY[2] - binningY[1]) / binningY[0] ):
        kwargs = {}
        kwargs[args.variables[0]] = var1
        kwargs[args.variables[1]] = var2

        # uncertainties
        c.reset()
        c.addUncertainty('lumi',        'lnN')
        c.addUncertainty('JEC',         'lnN')
        c.addUncertainty('fake',        'lnN')

        #processes = [ "signal", "ttX_SM_yield" ] + [s.name for s in bg]
        signal_rate                  = {}
        for i_region, region in enumerate(regions):
            signal_rate[region] = w.get_weight_yield( ttX_coeffList[region], **kwargs) - ttX_SM_rate[region] 

            bin_name = "Region_%i" % i_region
            nice_name = region.__str__()
            c.addBin(bin_name,['ttX_SM'] + [s.name for s in bg], nice_name)
            c.specifyObservation( bin_name, observation[region] )

            c.specifyFlatUncertainty( 'lumi', 1.05 )

            c.specifyExpectation( bin_name, 'signal', signal_rate[region] )
            c.specifyUncertainty( 'JEC', bin_name, 'signal', signal_jec_uncertainty[region])
            c.specifyUncertainty( 'fake',bin_name, 'signal', signal_fakerate_uncertainty[region])

            c.specifyExpectation( bin_name, 'ttX_SM', ttX_SM_rate[region] )
            c.specifyUncertainty( 'JEC', bin_name, 'ttX_SM', ttX_SM_jec_uncertainty[region])
            c.specifyUncertainty( 'fake',bin_name, 'ttX_SM', ttX_SM_fakerate_uncertainty[region])

            for background in bg:
                c.specifyExpectation( bin_name, background.name, background_rate[region][background.name] )
                c.specifyUncertainty( 'JEC', bin_name, background.name, background_jec_uncertainty[region][background.name])
                c.specifyUncertainty( 'fake',bin_name, background.name, background_fakerate_uncertainty[region][background.name])
                
        c.writeToFile( './tmp_limit_card.txt' ) 

        # try to adjust rmax with some margin
        exp_tot_sigmas = 0
        max_rmax = float('inf')
        for region in regions:
            tot_background = sum( [ background_rate[region][background.name] for background in bg ] )
            exp_tot_sigmas += abs(signal_rate[region])/sqrt( tot_background )
            # avoid total neg. yield
            if signal_rate[region]<0: 
                max_r = -tot_background/signal_rate[region]
                if max_r<max_rmax:
                    max_rmax = max_r
            
        rmax_est = 400./exp_tot_sigmas if exp_tot_sigmas>0 else 1.
        if max_rmax < rmax_est:
            rmax_est = 0.9*max_rmax # safety margin such that at least +10% total yield survives in the smallest SR

        profiledLoglikelihoodFit = ProfiledLoglikelihoodFit( './tmp_limit_card.txt' )
        profiledLoglikelihoodFit.make_workspace(rmin=0, rmax=rmax_est)
        #expected_limit = profiledLoglikelihoodFit.calculate_limit( calculator = "frequentist" )
        expected_limit = profiledLoglikelihoodFit.calculate_limit( calculator = "asymptotic" )
        logger.info( "Expected Limit: %f", expected_limit[0] )
        limit.SetBinContent( limit.FindBin(var1, var2), expected_limit[0] )
        profiledLoglikelihoodFit.cleanup()

limitPlot = Plot2D.fromHisto( '_'.join(args.variables), texX = args.variables[0], texY = args.variables[1], histos = [[limit]])
limitPlot.drawOption = "colz"
ROOT.gStyle.SetPaintTextFormat("2.2f")

WC_directory = '_'.join(args.parameters).rstrip('0').replace('-','m').replace('.','p') if args.parameters is not None else 'SM'
plot_directory_ = os.path.join(\
    plot_directory,
    '%s_%s'%(args.level, args.version),
    signal.name,
    'limit_small' if args.small else 'limit',
    WC_directory)

plotting.draw2D( limitPlot, plot_directory = plot_directory_, extensions = ["png", "root", "pdf"])

