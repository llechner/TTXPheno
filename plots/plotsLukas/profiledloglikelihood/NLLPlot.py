''' Plot script WC parameter LogLikelihood
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
argParser.add_argument('--version',         action='store',     default='v7', help='Appendix to plot directory')
argParser.add_argument('--process',         action='store',     default='ttZ_3l', nargs='?', choices=['ttZ_3l', 'ttZ_4l', 'ttgamma_1l', 'ttgamma_2l'], help="which process to calculate?")
argParser.add_argument('--sample',          action='store',     default='fwlite_ttZ_ll_LO_order2_15weights_ref', help='Sample name specified in sample/python/benchmarks.py, e.g. fwlite_ttZ_ll_LO_order2_15weights_ref')
argParser.add_argument('--order',           action='store',     default=2, help='Polynomial order of weight string (e.g. 2)')
argParser.add_argument('--selection',       action='store',     default='lepSel3-onZ-njet3p-nbjet1p-Zpt0', help="Specify cut.")
argParser.add_argument('--small',           action='store_true', help='Run only on a small subset of the data?')
argParser.add_argument('--level',           action='store',     default='gen', nargs='?', choices=['reco', 'gen'], help='Which level of reconstruction? reco, gen')
argParser.add_argument('--variables' ,      action='store',     default = ['cpQM', 'cpt'], type=str, nargs='+', help = "argument plotting variables")
argParser.add_argument('--binning',         action='store',     default = [12, -8, 40, 12, -30, 18], type=int, nargs='+', help = "argument parameters")
argParser.add_argument('--luminosity',      action='store',     default=150, help='Luminosity for weighting the plots')

args = argParser.parse_args()

if len(args.binning) != 6:
    raise ValueError('Binning must be 6 values in the form BINSX STARTX STOPX BINSY STARTY STOPY! Input: %s' %' '.join(args.binning))

binningX = args.binning[:3]
binningY = args.binning[3:]

if (binningX[2] - binningX[1]) % binningX[0] != 0:
    raise ValueError('Binning: Difference of upper and lower bound must be a multiple of the number of bins: (%s - %s) / %s != integer'%(str(binningX[2]), str(binningX[1]), str(binningX[0])) )

if (binningY[2] - binningY[1]) % binningY[0] != 0:
    raise ValueError('Binning: Difference of upper and lower bound must be a multiple of the number of bins: (%s - %s) / %s != integer'%(str(binningY[2]), str(binningY[1]), str(binningY[0])) )

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
sample_file     = "$CMSSW_BASE/python/TTXPheno/samples/benchmarks.py"
loadedSamples   = imp.load_source( "samples", os.path.expandvars( sample_file ) )

ttXSample       = getattr( loadedSamples, args.sample )
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

if args.process == 'ttZ_3l': bg = [ WZSample, tWZSample, tZqSample, ttgammaSample ]
elif args.process == 'ttZ_4l': bg = [ WZSample, tWZSample, tZqSample, ttgammaSample ]
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

if args.small:
    for s in [ttXSample] + bg:
        s.reduceFiles( to = 5 )

# configure samples
for s in [ttXSample] + bg:

    s.event_factor = s.nEvents / float( s.chain.GetEntries() )
    s.weightInfo = WeightInfo( s.reweight_pkl )
    s.weightInfo.set_order( args.order )
    s.setSelectionString( selectionString )

    if checkReferencePoint( s ):
#        s.read_variables = ["ref_lumiweight1fb/F", VectorTreeVariable.fromString('p[C/F]', nMax=2000) ]
        s.setWeightString( 'ref_lumiweight1fb*(%s)*(%s)'%( str(args.luminosity), str(s.event_factor) ) )
    else:
#        s.read_variables = ["lumiweight1fb/F"]
        s.setWeightString( 'lumiweight1fb*(%s)*(%s)'%( str(args.luminosity), str(s.event_factor) ) )

catPhoton_variables = [ "signalPhoton/I", "isrPhoton/I", "lepPhoton/I", "nonIsoPhoton/I", "fakePhoton/I"]
# overlap removal and signal selection
if args.process.split('_')[0] == 'ttgamma':
#    ttXSample.read_variables        += catPhoton_variables
#    ttgammaIsrSample.read_variables += catPhoton_variables
#    ttSample.read_variables         += catPhoton_variables

    ttXSample.addSelectionString(        "signalPhoton==1"                 ) #cat a1 <- the one and only signal
    ttgammaIsrSample.addSelectionString( "isrPhoton==1"                    ) #cat a2
    ttSample.addSelectionString(         "(signalPhoton!=1&&isrPhoton!=1)" ) #cat b,c,d

for var in args.variables:
    if var not in ttXSample.weightInfo.variables:
        raise ValueError('Input variable not in gridpack: %s' %var)

observation                  = {}

signal_jec_uncertainty       = {}
signal_fakerate_uncertainty  = {}

ttX_SM_rate                  = {}
#ttX_SM_jec_uncertainty       = {}
#ttX_SM_fakerate_uncertainty  = {}
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

#    ttX_SM_jec_uncertainty      [region] = 1.05
#    ttX_SM_fakerate_uncertainty [region] = 1.0  # signal has no FR uncertainty

    # signal uncertainties
    signal_jec_uncertainty      [region] = 1.05
    signal_fakerate_uncertainty [region] = 1.0  # signal has no FR uncertainty

    background_rate[region]                 = {}
    background_fakerate_uncertainty[region] = {}
    background_jec_uncertainty[region]      = {}

    for i_background, background in enumerate(bg):
        # compute bg yield for this region (this is the final code)

        background.setSelectionString( '&&'.join( [ selectionString, region.cutString() ] ) )

        background_rate                 [region][background.name] = background.getYieldFromDraw()['val']
        background_fakerate_uncertainty [region][background.name] = 1   + 0.03*i_region*(i_background+1) #change that
        background_jec_uncertainty      [region][background.name] = 1.2 - 0.02*i_region*(i_background+1) #change that

    # Our expected observation :-)
    observation[region] = int( sum( background_rate[region].values() ) + ttX_SM_rate[region] )


# Write temporary card file
from TTXPheno.Tools.cardFileWriter import cardFileWriter
c = cardFileWriter.cardFileWriter()

# Limit plot
from TTXPheno.Analysis.ProfiledLoglikelihoodFit import ProfiledLoglikelihoodFit

nll_plot = ROOT.TH2F( 'nll_plot', 'nll_plot', *(binningX + binningY) )

for var1 in range( binningX[1], binningX[2], ( binningX[2] - binningX[1]) / binningX[0] ):
    for var2 in range( binningY[1], binningY[2], ( binningY[2] - binningY[1]) / binningY[0] ):

        kwargs = { args.variables[0]:var1, args.variables[1]:var2 }

        # uncertainties
        c.reset()
        c.addUncertainty('lumi',        'lnN')
        c.addUncertainty('JEC',         'lnN')
        c.addUncertainty('fake',        'lnN')

        signal_rate                  = {}
        for i_region, region in enumerate(regions):

            signal_rate[region] = ttXSample.weightInfo.get_weight_yield( ttX_coeffList[region], **kwargs)

            bin_name = "Region_%i" % i_region
            nice_name = region.__str__()
            c.addBin(bin_name, ['_'.join(s.name.split('_')[1:3]) for s in bg], nice_name)
            c.specifyObservation( bin_name, observation[region] )

            c.specifyFlatUncertainty( 'lumi', 1.05 )

            c.specifyExpectation( bin_name, 'signal', signal_rate[region] )
            c.specifyUncertainty( 'JEC', bin_name, 'signal', signal_jec_uncertainty[region])
            c.specifyUncertainty( 'fake',bin_name, 'signal', signal_fakerate_uncertainty[region])

            #c.specifyExpectation( bin_name, 'ttX_SM', ttX_SM_rate[region] )
            #c.specifyUncertainty( 'JEC', bin_name, 'ttX_SM', ttX_SM_jec_uncertainty[region])
            #c.specifyUncertainty( 'fake',bin_name, 'ttX_SM', ttX_SM_fakerate_uncertainty[region])

            for background in bg:
                c.specifyExpectation( bin_name, '_'.join( background.name.split('_')[1:3] ), background_rate[region][background.name] )
                c.specifyUncertainty( 'JEC', bin_name, '_'.join( background.name.split('_')[1:3] ), background_jec_uncertainty[region][background.name])
                c.specifyUncertainty( 'fake',bin_name, '_'.join( background.name.split('_')[1:3] ), background_fakerate_uncertainty[region][background.name])
                
        nameList = ttXSample.name.split('_')[1:3] + args.variables + args.binning + [ args.level, args.version, args.order, args.luminosity, args.selection, 'small' if args.small else 'full' ]
        cardname = '%s_nll_card'%'_'.join( map( str, nameList ) )
        c.writeToFile( './tmp/%s.txt'%cardname )

        profiledLoglikelihoodFit = ProfiledLoglikelihoodFit( './tmp/%s.txt'%cardname )
        profiledLoglikelihoodFit.make_workspace(rmin=0, rmax=1)
        #expected_limit = profiledLoglikelihoodFit.calculate_limit( calculator = "frequentist" )
        nll = profiledLoglikelihoodFit.likelihoodTest()
        logger.info( "NLL: %f", nll)
        nll_plot.SetBinContent( nll_plot.FindBin(var1, var2), nll )
        profiledLoglikelihoodFit.cleanup()

nll_SM = nll_plot.GetBinContent( nll_plot.FindBin(0,0) )
for bin_x in range(nll_plot.GetNbinsX()+1):
    for bin_y in range(nll_plot.GetNbinsY()+1):
        nll_plot.SetBinContent(bin_x, bin_y, nll_plot.GetBinContent(  nll_plot.GetBin(bin_x, bin_y) ) - nll_SM)

nllPlot = Plot2D.fromHisto('_'.join(args.variables), texX = args.variables[0], texY = args.variables[1], histos = [[nll_plot]])
nllPlot.drawOption = "colz"
ROOT.gStyle.SetPaintTextFormat("2.2f")        

plot_directory_ = os.path.join(\
    plot_directory,
    '%s_%s'%(args.level, args.version),
    ttXSample.name,
    'NLL_small' if args.small else 'NLL')

plotting.draw2D( nllPlot, plot_directory = plot_directory_, extensions = ["png"])

