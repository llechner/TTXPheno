''' Fit EFT parameters to observation.
    Currently, this class is a placeholder to interface to Wolfgangs code.
'''

# Standard imports 
import sys
import ROOT
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
from TTXPheno.Tools.cutInterpreterGen import cutInterpreter
from TTXPheno.samples.benchmarks import * 
from TTXPheno.Tools.user import plot_directory

# Sample
ttZ_sample = fwlite_ttZ_ll_LO_order2_15weights_ref

# reduce dataset
ttZ_sample.reduceFiles( to = 5 )
# approximately compensate reduction of files
ttZ_sample.setWeightString("40")

# get the reweighting function
from TTXPheno.Tools.WeightInfo import WeightInfo
ttZ_sample.weightInfo = WeightInfo(ttZ_sample.reweight_pkl)
ttZ_sample.weightInfo.set_order( 2 )

# Load the analysis regions
from TTXPheno.Analysis.regions import regions

class FakeBackground:
    ''' Just a fake background number producer. Ignores everything, returns a number.
    '''
    def __init__( self, name, relative_fraction):
        self.name              = name
        self.relative_fraction = relative_fraction
    def setSelectionString( self, *args, **kwargs):
        pass 

# Lumi
lumi        =  150

# Backgrounds
backgrounds = [FakeBackground("bkg1",.1), FakeBackground("bkg2", .2)]

# set selection string
selectionString = cutInterpreter.cutString('lepSel3-onZ-njet3p-nbjet1p')
for sample in [ttZ_sample]+list(backgrounds):
    sample.setSelectionString( cutInterpreter.cutString('lepSel3-onZ-njet3p-nbjet1p') )

# parameter point
params = {'ctZ': 0., 'ctZI':0.}     # Can specify any EFT parameter here

observation                  = {}
signal_jec_uncertainty       = {}
signal_fakerate_uncertainty  = {}
ttZ_SM_rate                  = {}
#ttZ_SM_jec_uncertainty       = {}
#ttZ_SM_fakerate_uncertainty  = {}
background_rate                 = {}
background_jec_uncertainty      = {}
background_fakerate_uncertainty = {}

ttZ_coeffList                = {}
for i_region, region in enumerate(regions):

    logger.info( "At region %s", region )

    # compute signal yield for this region (this is the final code)
    ttZ_coeffList[region] = ttZ_sample.weightInfo.getCoeffListFromDraw( ttZ_sample, region.cutString(), weightString='150*ref_lumiweight1fb' ) 
    ## TTZ SM
    ttZ_SM_rate[region] = ttZ_sample.weightInfo.get_weight_yield( ttZ_coeffList[region] )
    #ttZ_SM_jec_uncertainty      [region] = 1.05 
    #ttZ_SM_fakerate_uncertainty [region] = 1.0  # signal has no FR uncertainty

    # signal uncertainties
    signal_jec_uncertainty      [region] = 1.05 
    signal_fakerate_uncertainty [region] = 1.0  # signal has no FR uncertainty

    background_rate[region]                = {}
    background_fakerate_uncertainty[region] = {}
    background_jec_uncertainty[region]      = {}
    for i_background, background in enumerate(backgrounds):
        # compute some  background yield
        background_rate[region][background.name] = background.relative_fraction*ttZ_SM_rate[region]

        # Guess the uncertainty
        background_fakerate_uncertainty[region][background.name] =  1   + 0.03*i_region*(i_background+1) 
        background_jec_uncertainty[region][background.name]      =  1.2 - 0.02*i_region*(i_background+1) 

    # Our expected observation :-)
    observation[region] = int( sum( background_rate[region].values() ) + ttZ_SM_rate[region] )


# Write temporary card file
from TTXPheno.Tools.cardFileWriter import cardFileWriter
c = cardFileWriter.cardFileWriter()

# Limit plot
from TTXPheno.Analysis.ProfiledLoglikelihoodFit import ProfiledLoglikelihoodFit

#binningX = binningY = [12, -12, 12]

binningX = [24, -24-6, 24-6]
binningY = [24, -24+16, 24+16]

nll_plot = ROOT.TH2F( 'nll_plot', 'nll_plot', *(binningX + binningY) )

for cpt in range( binningX[1], binningX[2], ( binningX[2] - binningX[1]) / binningX[0] ):
    for cpQM in range( binningY[1], binningY[2], ( binningY[2] - binningY[1]) / binningY[0] ):
#for cpt in [4]:
#    for cpQM in [4]:
        kwargs = {'cpt':cpt, 'cpQM':cpQM}

        # uncertainties
        c.reset()
        c.addUncertainty('lumi',        'lnN')
        c.addUncertainty('JEC',         'lnN')
        c.addUncertainty('fake',        'lnN')

        signal_rate                  = {}
        for i_region, region in enumerate(regions):
            signal_rate[region] = ttZ_sample.weightInfo.get_weight_yield( ttZ_coeffList[region], **kwargs) 

            bin_name = "Region_%i" % i_region
            nice_name = region.__str__()
            c.addBin(bin_name,[s.name for s in backgrounds], nice_name)
            c.specifyObservation( bin_name, observation[region] )

            c.specifyFlatUncertainty( 'lumi', 1.05 )

            c.specifyExpectation( bin_name, 'signal', signal_rate[region] )
            c.specifyUncertainty( 'JEC', bin_name, 'signal', signal_jec_uncertainty[region])
            c.specifyUncertainty( 'fake',bin_name, 'signal', signal_fakerate_uncertainty[region])

            #c.specifyExpectation( bin_name, 'ttZ_SM', ttZ_SM_rate[region] )
            #c.specifyUncertainty( 'JEC', bin_name, 'ttZ_SM', ttZ_SM_jec_uncertainty[region])
            #c.specifyUncertainty( 'fake',bin_name, 'ttZ_SM', ttZ_SM_fakerate_uncertainty[region])

            for background in backgrounds:
                c.specifyExpectation( bin_name, background.name, background_rate[region][background.name] )
                c.specifyUncertainty( 'JEC', bin_name, background.name, background_jec_uncertainty[region][background.name])
                c.specifyUncertainty( 'fake',bin_name, background.name, background_fakerate_uncertainty[region][background.name])
                
        c.writeToFile( './tmp_nll_card.txt' ) 

        profiledLoglikelihoodFit = ProfiledLoglikelihoodFit( './tmp_nll_card.txt' )
        profiledLoglikelihoodFit.make_workspace(rmin=0, rmax=1)
        #expected_limit = profiledLoglikelihoodFit.calculate_limit( calculator = "frequentist" )
        nll = profiledLoglikelihoodFit.likelihoodTest()
        logger.info( "NLL: %f", nll)
        nll_plot.SetBinContent( nll_plot.FindBin(cpt, cpQM), nll )
        profiledLoglikelihoodFit.cleanup()

nll_SM = nll_plot.GetBinContent( nll_plot.FindBin(0,0))
for bin_x in range(nll_plot.GetNbinsX()+1):
    for bin_y in range(nll_plot.GetNbinsY()+1):
        nll_plot.SetBinContent(bin_x, bin_y, nll_plot.GetBinContent(  nll_plot.GetBin(bin_x, bin_y) ) - nll_SM)

nllPlot = Plot2D.fromHisto( "2D_nll_2", texX = "cpt", texY = "cpQM", histos = [[nll_plot]])
nllPlot.drawOption = "colz"
ROOT.gStyle.SetPaintTextFormat("2.2f")        
plotting.draw2D( nllPlot, plot_directory = os.path.join( plot_directory, 'nll', ttZ_sample.name), extensions = ["png"])

