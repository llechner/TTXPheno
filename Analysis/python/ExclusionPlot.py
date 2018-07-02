''' Fit EFT parameters to observation.
    Currently, this class is a placeholder to interface to Wolfgangs code.
'''

class FakeBackground:
    ''' Just a fake background number producer. Ignores everything, returns a number.
    '''
    def __init__( self, name, relative_fraction):
        self.name              = name
        self.relative_fraction = relative_fraction
    def setSelectionString( self, *args, **kwargs):
        pass 


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

# Sample
ttZ_sample = fwlite_ttZ_ll_LO_order2_15weights_ref

# reduce dataset
ttZ_sample.reduceFiles( to = 1 )
# approximately compensate reduction of files
ttZ_sample.setWeightString("200")

# get the reweighting function
from TTXPheno.Tools.WeightInfo import WeightInfo, getCoeffListFromDraw
ttZ_sample.weightInfo = WeightInfo(ttZ_sample.reweight_pkl)
ttZ_sample.weightInfo.set_order( 2 )

# Load the analysis regions
from TTXPheno.Analysis.regions import regions

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

ttZ_SM_yield                 = {}
signal_yield                 = {}
jec_uncertainty_signal       = {}
fakerate_uncertainty_signal  = {}
rates                        = {}
observation                  = {}
jec_uncertainty              = {}
fakerate_uncertainty         = {}
ttZ_coeffList                = {}
for i_region, region in enumerate(regions):

    logger.info( "At region %s", region )

    # compute signal yield for this region (this is the final code)
    ttZ_coeffList[region] = getCoeffListFromDraw( ttZ_sample, region.cutString(), weightString='150*ref_lumiweight1fb' ) 
    # TTZ SM
    ttZ_SM_rate[region] = get_weight_yield( ttZ_coeffList[region] )
    ttZ_SM_jec_uncertainty      [region] = 1.05 
    ttZ_SM_fakerate_uncertainty [region] = 1.0  # signal has no FR uncertainty

    # temporary, replace later
    signal_rate                 [region] = 0. 
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
    observation[region] = int( sum( background_rate[region] ) + ttZ_SM_rate[region] )


# Write temporary card file
from TopEFT.Tools.cardFileWriter import cardFileWriter
c = cardFileWriter.cardFileWriter()
cardFileName = os.path.join( 'tmp_card.txt')
# uncertainties
c.reset()
c.addUncertainty('lumi',        'lnN')
c.addUncertainty('JEC',         'lnN')
c.addUncertainty('fake',        'lnN')

#processes = [ "signal", "ttZ_SM_yield" ] + [s.name for s in backgrounds]

for i_region, region in enumerate(regions):
    bin_name = "Region_%i" % i_region
    nice_name = region.__str__()
    c.addBin(bin_name,['signal', 'ttZ_SM'] + [s.name for s in backgrounds], nice_name)
    c.specifyObservation( bin_name, observation[region] )

    c.specifyFlatUncertainty( 'lumi', 1.05 )

    c.specifyExpectation(binname, 'signal', signal_rate['region'] )
    c.specifyUncertainty( 'JEC', bin_name, 'signal', signal_jec_uncertainty[region])
    c.specifyUncertainty( 'fake',bin_name, 'signal', signal_fakerate_uncertainty[region])

    c.specifyExpectation(binname, 'ttZ_SM', ttZ_SM_rate['region'] )
    c.specifyUncertainty( 'JEC', bin_name, 'ttZ_SM', ttZ_SM_jec_uncertainty[region])
    c.specifyUncertainty( 'fake',bin_name, 'ttZ_SM', ttZ_SM_fakerate_uncertainty[region])

    for background in backgrounds:
        c.specifyExpectation(binname, background.name, background_rate['region'] )
        c.specifyUncertainty( 'JEC', bin_name, background.name, background_jec_uncertainty[region])
        c.specifyUncertainty( 'fake',bin_name, background.name, background_fakerate_uncertainty[region])
        

