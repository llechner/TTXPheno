''' Fit EFT parameters to observation.
    Currently, this class is a placeholder to interface to Wolfgangs code.
'''

class LikelihoodFit:

    def __init__(self, lumi, backgrounds, signal, selectionString, regions):

        self.lumi        = lumi
        self.backgrounds = backgrounds
        self.signal      = signal
    
        # Apply selection
        for background in self.backgrounds:
            background.setSelectionString( selectionString )

        self.regions = regions

    def run( self):
   
        signal_yield          = {}
        background_yields     = {}
        observation           = {}
        jec_uncertainty       = {}
        fakerate_uncertainty  = {}
 
        for i_region, region in enumerate(regions):

            logger.info( "At region %s", region )

            # parameter point
            params = {'ctZ': 0., 'ctZI':0.}     # Can specify any EFT parameter here

            # compute signal yield for this region (this is the final code)
            signal_yield[region] = self.signal.getYieldFromDraw( 
                  region.cutString(), 
                  weightString = "%f*ref_lumiweight1fb*(%s)"%(self.lumi, signal.weightInfo.get_weight_string(**params))
                )['val']

            # compute some  background yield
            background_yields[region] = [ background.relative_fraction*signal_yield[region] for background in self.backgrounds ]

            # Our expected observation :-)
            observation[region] = int( sum( background_yields[region] ) + signal_yield[region] )
    
            # Guess the uncertainty
            fakerate_uncertainty[region] = 1 + 0.1*i_region
            jec_uncertainty[region] = 1.3 - 0.05*i_region

        # Wolfgang, here comes your data!
        logger.info("Now I perform the fit! My data is:")
        for region in regions:
            logger.info( "Obs: %5i Signal: %5.2f Bkg yields: %s jec: %5.2f fake: %5.2f", 
                observation[region], 
                signal_yield[region],
                ",".join( map("{:7.2f}".format,  background_yields[region] ) ),
                jec_uncertainty[region],
                fakerate_uncertainty[region]
            )

class FakeBackground:
    ''' Just a fake background number producer. Ignores everything, returns a number.
    '''
    def __init__( self, relative_fraction):
        self.relative_fraction = relative_fraction

    def setSelectionString( self, *args, **kwargs):
        pass 


if __name__ == '__main__':
    
    # RootTools
    from RootTools.core.standard import *

    # Logger
    import TTXPheno.Tools.logger as logger
    import RootTools.core.logger as logger_rt
    logger    = logger.get_logger(   'DEBUG', logFile = None)
    logger_rt = logger_rt.get_logger('INFO', logFile = None)

    # TTXPheno
    from TTXPheno.Tools.cutInterpreter import cutInterpreter
    from TTXPheno.samples.benchmarks import * 

    # Sample
    signal = fwlite_ttZ_ll_LO_order2_15weights_ref

    # reduce dataset
    signal.reduceFiles( to = 1 )
    # approximately compensate reduction of files
    signal.setWeightString("200")

    # get the reweighting function
    from TTXPheno.Tools.WeightInfo import WeightInfo
    signal.weightInfo = WeightInfo(signal.reweight_pkl)
    signal.weightInfo.set_order( 2 )


    # Load the analysis regions
    from TTXPheno.Analysis.regions import regions

    # set selection string
    signal.setSelectionString( cutInterpreter.cutString('lepSel3-onZ-njet3p-nbjet1p') )

    fit = LikelihoodFit(
            150,                                                        # Lumi 
            [FakeBackground(.1), FakeBackground(.2)],                   # background samples -> 10% and 20% of the signal
            signal,                                                     # signal
            cutInterpreter.cutString('lepSel3-onZ-njet3p-nbjet1p'),     # preselection
            regions                                                     # signal regions
        )

    fit.run()
