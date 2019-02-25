''' Plot script WC parameter LogLikelihood
'''

# Standard imports 
import sys
import ROOT
import imp
import pickle
import ctypes
import numpy as np
import itertools
import operator

from math import sqrt
# turn off graphics
ROOT.gROOT.SetBatch( True )

# RootTools
from RootTools.core.standard import *

from plot_helpers import getUncertaintyValue, getObservationValue
from multiprocessing import Pool

# Logger
import TTXPheno.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   'DEBUG', logFile = None)
logger_rt = logger_rt.get_logger('INFO', logFile = None)

# TTXPheno
from TTXPheno.samples.benchmarks import * 
from TTXPheno.Tools.user import plot_directory, cardfileLocation

# get the reweighting function
from TTXPheno.Tools.WeightInfo import WeightInfo

ROOT.gStyle.SetNumberContours(255)

# Arguments
import argparse

argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--version',            action='store',     default='test', help='Appendix to plot directory')
argParser.add_argument('--process',            action='store',     default='ttZ_3l', nargs='?', choices=['ttZ_3l', 'ttZ_4l', 'ttgamma_1l', 'ttgamma_2l'], help="which process to calculate?")
argParser.add_argument('--bestFit',            action='store_true', help='Run combine with bestFit scenario (wide r ranges)')
argParser.add_argument('--removeCardFiles',    action='store_true', help='remove cardfiles after calculation?')
argParser.add_argument('--useCombine',         action='store_true', help='use Higgs Combine Tool?')
argParser.add_argument('--fitOnly',            action='store_true', help='do you already have the cardfiles?')
argParser.add_argument('--sample',             action='store',     default='fwlite_ttZ_ll_LO_order2_15weights_ref', help='Sample name specified in sample/python/benchmarks.py, e.g. fwlite_ttZ_ll_LO_order2_15weights_ref')
argParser.add_argument('--order',              action='store',     default=2, help='Polynomial order of weight string (e.g. 2)')
argParser.add_argument('--selection',          action='store',     default='lepSel3-onZ-njet3p-nbjet1p-Zpt0', help="Specify cut.")
argParser.add_argument('--small',              action='store_true', help='Run only on a small subset of the data?')
argParser.add_argument('--contours',           action='store_true', help='draw 1sigma and 2sigma contour line?')
argParser.add_argument('--smooth',             action='store_true', help='smooth histogram?')
argParser.add_argument('--level',              action='store',     default='reco', nargs='?', choices=['reco', 'gen'], help='Which level of reconstruction? reco, gen')
argParser.add_argument('--variable' ,         action='store',     default = 'ctZ', type=str, nargs='?', help = "argument plotting variable")
argParser.add_argument('--binning',            action='store',     default = [1, -2, 2], type=float, nargs=3, help = "argument parameters")
argParser.add_argument('--yRange',             action='store',     default = [None, None], type=float, nargs=2, help = "argument parameters")
argParser.add_argument('--luminosity',         action='store',     default=150, type=int, help='Luminosity for weighting the plots')
argParser.add_argument('--cores',              action='store',     default=8, type=int, help='number of cpu cores for multicore processing')
argParser.add_argument('--overwrite',          action='store_true', help='overwrite datafile?')
argParser.add_argument('--binMultiplier',      action='store',     default=3, type=int, help='bin multiplication factor')
argParser.add_argument('--detector',           action='store',     default='CMS', nargs='?', choices=['CMS', 'ATLAS', 'phase2_CMS'], help='Which Delphes detector simulation?')
argParser.add_argument('--scale14TeV',         action='store_true', help='scale 13 TeV cross-sections to 14 Tev?')
argParser.add_argument('--additionalCardFile', action='store',     default='TopEFTCardFile.txt', help='Cardfile where additional uncertainties are taken from')
argParser.add_argument('--addNonPrompt',       action='store_true', help='add nonPrompt?')
argParser.add_argument('--addUncertainties',   action='store',     default = ['trigger_2016','scale','scale_sig','PDF','PartonShower','nonprompt','WZ_xsec', 'ttX'], type=str, help = "add additional uncertainties from cardFile")
argParser.add_argument('--addBinNumberShift',  action='store',     default = 0, type=int, help = "which bin number does the region start in the additional card file?")
argParser.add_argument('--uncertaintyScale',   action='store',     default = 0.5, type=float, help = "scale factor for additional uncertainties")
argParser.add_argument('--statOnly',           action='store_true', help='use only statistical uncertainties')
argParser.add_argument('--noExpUnc',          action='store_true', help='use only statistical and theory uncertainties')

args = argParser.parse_args()

for unc in ['trigger_2016', 'nonprompt']:
    if unc in args.addUncertainties and args.noExpUnc: args.addUncertainties.remove(unc)

if args.level == 'gen':
    if 'ttZ' in args.process.split('_'):
        if args.small: from TTXPheno.Analysis.regions import genttZRegionsSmall as regions
        else:          from TTXPheno.Analysis.regions import genttZRegions as regions
    elif 'ttgamma' in args.process.split('_'):
        if args.small: from TTXPheno.Analysis.regions import genttgammaRegionsSmall as regions
        else:          from TTXPheno.Analysis.regions import genttgammaRegions as regions

    # Import additional functions/classes specified for the level of reconstruction
    from TTXPheno.Tools.cutInterpreterGen import cutInterpreter

elif args.level == 'reco':
    if 'ttZ' in args.process.split('_'):
        if args.small: from TTXPheno.Analysis.regions import recottZRegionsSmall as regions
        else:          from TTXPheno.Analysis.regions import recottZRegions as regions
    elif 'ttgamma' in args.process.split('_'):
        if args.small: from TTXPheno.Analysis.regions import recottgammaRegionsSmall as regions
        else:          from TTXPheno.Analysis.regions import recottgammaRegions as regions

    # Import additional functions/classes specified for the level of reconstruction
    from TTXPheno.Tools.cutInterpreterReco import cutInterpreter

if args.binning[0] > 1:
    xRange = np.linspace( args.binning[1], args.binning[2], int(args.binning[0]), endpoint=False)
    xRange = [ el + 0.5 * ( xRange[1] - xRange[0] ) for el in xRange ]
else:
    xRange = [ 0.5 * ( args.binning[1] + args.binning[2] ) ]

addon = []
if args.statOnly: addon += ["statOnly"]
if args.noExpUnc: addon += ["noExpUnc"]
#save data file
filename = '_'.join( ['nll', args.detector ] + args.sample.split('_')[1:3] + [args.variable] + map( str, args.binning ) + [ args.selection, str(args.luminosity), "14TeV" if args.scale14TeV else "13TeV" ] + addon ) + '.data'

#do the calculation
if not os.path.isfile('dat/' + filename) or args.overwrite:

    if not args.fitOnly:

        # Import samples
        sample_file     = "$CMSSW_BASE/python/TTXPheno/samples/benchmarks.py"
        loadedSamples   = imp.load_source( "samples", os.path.expandvars( sample_file ) )

        ttXSample       = getattr( loadedSamples, args.sample + '_%s' %args.detector )
        WZSample        = getattr( loadedSamples, 'fwlite_WZ_lep_LO_order2_15weights_%s' %args.detector )
        tWZSample       = getattr( loadedSamples, 'fwlite_tWZ_LO_order2_15weights_%s' %args.detector )
        tZqSample       = getattr( loadedSamples, 'fwlite_tZq_LO_order2_15weights_%s' %args.detector )
        ttgammaSample   = getattr( loadedSamples, 'fwlite_ttgamma_bg_LO_order2_15weights_%s' %args.detector )

        #if args.process.split('_')[0] == 'ttgamma':
        #    ttgammaIsrSample  = copy.deepcopy( ttXSample ) #select ttgamma events with isolated gamma from ISR (cat a2)
        #    ttgammaIsrSample.name = 'fwlite_ttgamma_ISR_LO_order2_15weights_ref'

        if args.process == 'ttZ_3l': bg = [ WZSample, tWZSample, tZqSample, ttgammaSample ]

        def checkReferencePoint( sample ):
            ''' check if sample is simulated with a reference point
            '''
            return pickle.load(file(sample.reweight_pkl))['ref_point'] != {}

        # set selection string
        selectionString      = cutInterpreter.cutString(args.selection)
        selectionString_up   = selectionString.replace('nBTag','nBTag_JEC_up').replace('nrecoJet','nrecoJets_JEC_up')
        selectionString_down = selectionString.replace('nBTag','nBTag_JEC_down').replace('nrecoJet','nrecoJets_JEC_down')

        # somehow has to be separate from the next loop
        if args.small:
            for s in [ttXSample] + bg:
                s.reduceFiles( to = 15 )

        # configure samples
        for s in [ttXSample] + bg:

            s.event_factor = s.nEvents / float( s.chain.GetEntries() )
            s.xsecScaleFactor = s.xsec14 / s.xsec if args.scale14TeV else 1.
            s.weightInfo = WeightInfo( s.reweight_pkl )
            s.weightInfo.set_order( args.order )
            s.setSelectionString( selectionString )

            if checkReferencePoint( s ):
                s.setWeightString( 'ref_lumiweight1fb*(%s)*(%s)*(%s)'%( str(args.luminosity), str(s.event_factor), str(s.xsecScaleFactor) ) )
            else:
                s.setWeightString( 'lumiweight1fb*(%s)*(%s)*(%s)'%( str(args.luminosity), str(s.event_factor), str(s.xsecScaleFactor) ) )

        # overlap removal
        if args.process.split('_')[0] == 'ttgamma':
            ttXSample.addSelectionString( "(nonIsoPhoton!=1)" ) 
            ttSample.addSelectionString(  "(nonIsoPhoton==1)" ) 

        if args.variable not in ttXSample.weightInfo.variables:
            raise ValueError('Input variable not in gridpack: %s' %args.variable)

        observation                    = {}

        signal_btagging_uncertainty    = {}
        signal_mistagging_uncertainty  = {}
        signal_jes_uncertainty         = {}
        signal_electronId_uncertainty  = {}
        signal_muonId_uncertainty      = {}

        ttX_SM_rate                  = {}
        ttX_coeffList                = {}

        ttX_coeffList_reweighted_btagging   = {}
        ttX_coeffList_reweighted_mistagging = {}
        ttX_coeffList_reweighted_muonId     = {}
        ttX_coeffList_reweighted_electronId = {}
        ttX_coeffList_reweighted_jes_up     = {}
        ttX_coeffList_reweighted_jes_down   = {}

        background_rate                   = {}
        background_btagging_uncertainty   = {}
        background_mistagging_uncertainty = {}
        background_jes_uncertainty        = {}
        background_electronId_uncertainty = {}
        background_muonId_uncertainty     = {}

        nonPromptObservation              = {}

        for i_region, region in enumerate(regions):
            # compute signal yield for this region (this is the final code)

            logger.info( "At region %s", region )

            # ttX SM
            ttX_coeffList[region]  = ttXSample.weightInfo.getCoeffListFromDraw( ttXSample, selectionString = region.cutString() )
            ttX_SM_rate[region]    = ttXSample.weightInfo.get_weight_yield( ttX_coeffList[region] )

            if not args.statOnly:
                # uncertainty coeffLists
                ttX_coeffList_reweighted_btagging[region]   = ttXSample.weightInfo.getCoeffListFromDraw( ttXSample, selectionString = region.cutString(), weightString="reweight_BTag_B" )
                ttX_coeffList_reweighted_mistagging[region] = ttXSample.weightInfo.getCoeffListFromDraw( ttXSample, selectionString = region.cutString(), weightString="reweight_BTag_L" )
                ttX_coeffList_reweighted_muonId[region]     = ttXSample.weightInfo.getCoeffListFromDraw( ttXSample, selectionString = region.cutString(), weightString="reweight_id_mu" )
                ttX_coeffList_reweighted_electronId[region] = ttXSample.weightInfo.getCoeffListFromDraw( ttXSample, selectionString = region.cutString(), weightString="reweight_id_ele" )

                ttXSample.setSelectionString( selectionString_up )
                ttX_coeffList_reweighted_jes_up[region]     = ttXSample.weightInfo.getCoeffListFromDraw( ttXSample, selectionString = region.cutString() )
    
                ttXSample.setSelectionString( selectionString_down )
                ttX_coeffList_reweighted_jes_down[region]   = ttXSample.weightInfo.getCoeffListFromDraw( ttXSample, selectionString = region.cutString() )

                # reset selectionstring
                ttXSample.setSelectionString( selectionString )

            background_rate[region]                   = {}
            background_btagging_uncertainty[region]   = {}
            background_mistagging_uncertainty[region] = {}
            background_jes_uncertainty[region]        = {}
            background_muonId_uncertainty[region]     = {}
            background_electronId_uncertainty[region] = {}

            for i_background, background in enumerate(bg):
                # compute bg yield for this region (this is the final code)

                background_rate                 [region][background.name] = background.getYieldFromDraw( selectionString=region.cutString() )['val']

                if not args.statOnly and not args.noExpUnc:
                    #calculate btagging uncert.
                    background_rate_reweighted                                = background.getYieldFromDraw( selectionString=region.cutString(), weightString="reweight_BTag_B" )['val']
                    background_btagging_uncertainty [region][background.name] = 1 + (( background_rate_reweighted - background_rate[region][background.name] ) / background_rate[region][background.name]) if background_rate[region][background.name] > 0 else 1.

                    #calculate mistagging uncert.
                    background_rate_reweighted                                = background.getYieldFromDraw( selectionString=region.cutString(), weightString="reweight_BTag_L" )['val']
                    background_mistagging_uncertainty [region][background.name] = 1 + (( background_rate_reweighted - background_rate[region][background.name] ) / background_rate[region][background.name]) if background_rate[region][background.name] > 0 else 1.

                    #calculate muon Id uncert.
                    background_rate_reweighted                                = background.getYieldFromDraw( selectionString=region.cutString(), weightString="reweight_id_mu" )['val']
                    background_muonId_uncertainty [region][background.name] = 1 + (( background_rate_reweighted - background_rate[region][background.name] ) / background_rate[region][background.name]) if background_rate[region][background.name] > 0 else 1.

                    background_rate_reweighted                                = background.getYieldFromDraw( selectionString=region.cutString(), weightString="reweight_id_ele" )['val']
                    #calculate electron Id uncert.
                    background_electronId_uncertainty [region][background.name] = 1 + (( background_rate_reweighted - background_rate[region][background.name] ) / background_rate[region][background.name]) if background_rate[region][background.name] > 0 else 1.

                    # set selectionstring to JES_up
                    background.setSelectionString( selectionString_up )
                    background_rate_reweighted_up                             = background.getYieldFromDraw( selectionString=region.cutString() )['val']
                    # set selectionstring to JES_up
                    background.setSelectionString( selectionString_down )
                    background_rate_reweighted_down                           = background.getYieldFromDraw( selectionString=region.cutString() )['val']
                    # reset selectionstring
                    background.setSelectionString( selectionString )
                    #calculate JES uncert.
                    background_jes_uncertainty      [region][background.name] = 1 + (( background_rate_reweighted_up - background_rate_reweighted_down ) / (2*background_rate[region][background.name])) if background_rate[region][background.name] > 0 else 1.

            nonPromptObservation[region] = 0.
            if args.addNonPrompt:
                # scale nonprompt observation value from Run2 to args.luminosity
                nonPromptObservation[region] = getObservationValue( args.additionalCardFile, args.addBinNumberShift + i_region, 'nonPromptDD' ) * float(args.luminosity) / 35.9

            # Our expected observation :-)
            # add nonPrompt observation to total observation
            observation[region] = int( round( sum( background_rate[region].values() ) + ttX_SM_rate[region] + nonPromptObservation[region] ) )


    # Write temporary card file
    from TTXPheno.Tools.cardFileWriter import cardFileWriter
#    c = cardFileWriter.cardFileWriter()
    if args.useCombine:
        from TTXPheno.Tools.user import combineReleaseLocation
#        c.releaseLocation = combineReleaseLocation
    else:
        # non CMS NLL plot
        from TTXPheno.Analysis.ProfiledLoglikelihoodFit import ProfiledLoglikelihoodFit

    def cuBWtoctWZ( cuB, cuW ):
        ''' transforms C_tZ and C_tW to C_uB and C_uW
            C_tZ = Re( -sW*C_uB + cW*C_uW )
            C_tW = Re( C_uW )
            arXiv: 1802.07237
        '''
        sW=0.4715
        cW=0.8819

        ctW = cuW
        ctZ = -sW*cuB + cW*cuW

        return ctZ, ctW


    def calculation( var ):
        kwargs = { args.variable:var }

        nameList = args.sample.split('_')[1:3] + [args.variable] + args.binning + [ args.level, args.version, args.order, args.luminosity, "14TeV" if args.scale14TeV else "13TeV", args.selection, 'small' if args.small else 'full', 'statOnly' if args.statOnly else 'fullUnc' if not args.noExpUnc else 'noExpUnc', var ]
        cardname = '%s_nll_card'%'_'.join( map( str, nameList ) )
        cardFilePath = os.path.join( cardfileLocation, cardname + '.txt' )

        c = cardFileWriter.cardFileWriter()
        if args.useCombine:
            c.releaseLocation = combineReleaseLocation

        if not args.fitOnly:
#            print 'run cardfile'

            # uncertainties
            c.reset()
            if not args.statOnly:
                if not args.noExpUnc:
                    c.addUncertainty('lumi',        'lnN')
                    c.addUncertainty('JES',         'lnN')
                    c.addUncertainty('btagging',    'lnN')
                    c.addUncertainty('mistagging',  'lnN')
                    c.addUncertainty('muonId',      'lnN')
                    c.addUncertainty('electronId',  'lnN')
                for unc in args.addUncertainties:
                    c.addUncertainty(unc,  'lnN')

            signal_rate                  = {}
            for i_region, region in enumerate(regions):

                signal_rate[region] = ttXSample.weightInfo.get_weight_yield( ttX_coeffList[region], **kwargs)

                if not args.statOnly and not args.noExpUnc:
                    # signal uncertainties
                    # btagging
                    signal_rate_reweighted   = ttXSample.weightInfo.get_weight_yield( ttX_coeffList_reweighted_btagging[region], **kwargs )
                    signal_btagging_uncertainty [region] = 1 + (( signal_rate_reweighted - signal_rate[region] ) / signal_rate[region]) if signal_rate[region] > 0 else 1.

                    # mistagging
                    signal_rate_reweighted   = ttXSample.weightInfo.get_weight_yield( ttX_coeffList_reweighted_mistagging[region], **kwargs )
                    signal_mistagging_uncertainty [region] = 1 + (( signal_rate_reweighted - signal_rate[region] ) / signal_rate[region]) if signal_rate[region] > 0 else 1.

                    # muonId
                    signal_rate_reweighted   = ttXSample.weightInfo.get_weight_yield( ttX_coeffList_reweighted_muonId[region], **kwargs )
                    signal_muonId_uncertainty [region] = 1 + (( signal_rate_reweighted - signal_rate[region] ) / signal_rate[region]) if signal_rate[region] > 0 else 1.

                    # electronId
                    signal_rate_reweighted   = ttXSample.weightInfo.get_weight_yield( ttX_coeffList_reweighted_electronId[region], **kwargs )
                    signal_electronId_uncertainty [region] = 1 + (( signal_rate_reweighted - signal_rate[region] ) / signal_rate[region]) if signal_rate[region] > 0 else 1.

                    # JES
                    signal_rate_reweighted_JES_up   = ttXSample.weightInfo.get_weight_yield( ttX_coeffList_reweighted_jes_up[region], **kwargs )
                    signal_rate_reweighted_JES_down = ttXSample.weightInfo.get_weight_yield( ttX_coeffList_reweighted_jes_down[region], **kwargs )
                    signal_jes_uncertainty[region] = 1 + (( signal_rate_reweighted_JES_up - signal_rate_reweighted_JES_down ) / (2*signal_rate[region])) if signal_rate[region] > 0 else 1.

                bin_name = "Region_%i" % i_region
                nice_name = region.__str__()
                c.addBin(bin_name, ['_'.join(s.name.split('_')[1:3]) for s in bg] + ['nonPrompt'] if args.addNonPrompt else ['_'.join(s.name.split('_')[1:3]) for s in bg], nice_name)

                c.specifyObservation( bin_name, observation[region] )

                c.specifyExpectation( bin_name, 'signal', signal_rate[region]                                 )

                if not args.statOnly:
                    if not args.noExpUnc:
                        c.specifyFlatUncertainty( 'lumi', 1.01 )
                        c.specifyUncertainty( 'JES',        bin_name, 'signal', signal_jes_uncertainty[region]        )
                        c.specifyUncertainty( 'btagging',   bin_name, 'signal', signal_btagging_uncertainty[region]   )
                        c.specifyUncertainty( 'mistagging', bin_name, 'signal', signal_mistagging_uncertainty[region] )
                        c.specifyUncertainty( 'muonId',     bin_name, 'signal', signal_muonId_uncertainty[region]     )
                        c.specifyUncertainty( 'electronId', bin_name, 'signal', signal_electronId_uncertainty[region] )

                    for unc in args.addUncertainties:
                        c.specifyUncertainty( unc,      bin_name, 'signal', 1+(getUncertaintyValue( args.additionalCardFile, args.addBinNumberShift + i_region, 'signal', unc )-1)*args.uncertaintyScale )

                if args.addNonPrompt:
                    # for nonpromt only nonpromt uncertainty is important
                    c.specifyExpectation( bin_name, 'nonPrompt', nonPromptObservation[region] )
                    if not args.statOnly: c.specifyUncertainty( 'nonprompt',      bin_name, 'nonPrompt', 1+(getUncertaintyValue( args.additionalCardFile, args.addBinNumberShift + i_region, 'nonPromptDD', 'nonprompt' )-1)*args.uncertaintyScale )

                #c.specifyExpectation( bin_name, 'ttX_SM', ttX_SM_rate[region] )
                #c.specifyUncertainty( 'JES', bin_name, 'ttX_SM', ttX_SM_jes_uncertainty[region])
                #c.specifyUncertainty( 'btagging',bin_name, 'ttX_SM', ttX_SM_btagging_uncertainty[region])

                for background in bg:
                    c.specifyExpectation( bin_name, '_'.join( background.name.split('_')[1:3] ), background_rate[region][background.name] )
                    if not args.statOnly:
                        if not args.noExpUnc:
                            c.specifyUncertainty( 'JES',        bin_name, '_'.join( background.name.split('_')[1:3] ), background_jes_uncertainty[region][background.name])
                            c.specifyUncertainty( 'btagging',   bin_name, '_'.join( background.name.split('_')[1:3] ), background_btagging_uncertainty[region][background.name])
                            c.specifyUncertainty( 'mistagging', bin_name, '_'.join( background.name.split('_')[1:3] ), background_mistagging_uncertainty[region][background.name])
                            c.specifyUncertainty( 'muonId',     bin_name, '_'.join( background.name.split('_')[1:3] ), background_muonId_uncertainty[region][background.name])
                            c.specifyUncertainty( 'electronId', bin_name, '_'.join( background.name.split('_')[1:3] ), background_electronId_uncertainty[region][background.name])
                        for unc in args.addUncertainties:
                            if 'tZq' in background.name.split('_') or 'ttgamma' in background.name.split('_') or 'tWZ' in background.name.split('_'): proc = 'TTX'
                            elif 'WZ' in background.name.split('_'): proc = 'WZ'
                            else: raise ValueError('Background not found: %s' %background.name)
                            c.specifyUncertainty( unc,      bin_name, '_'.join( background.name.split('_')[1:3] ), 1+(getUncertaintyValue( args.additionalCardFile, args.addBinNumberShift + i_region, proc, unc )-1)*args.uncertaintyScale )
                    
            c.writeToFile( cardFilePath )

        else:
            logger.info( "Running only NLL Fit with given CardFile %s"%cardFilePath)

        if not os.path.isfile( cardFilePath ):
            raise ValueError('CardFiles not found! Run script without --fitOnly!')

        if args.useCombine:
            # use the official cms combine tool
#                c.calcNuisances( cardFilePath, bestFit=args.bestFit )
            nll = c.calcNLL( cardFilePath, bestFit=args.bestFit )
#            nll = nll['nll0'] #pre-fit
            nll = nll['nll_abs'] #post-fit

            if args.removeCardFiles:
                for file in os.listdir( cardfileLocation ):
                    if file.startswith( cardname ):
                        os.remove( os.path.join( cardfileLocation, file ) )

        else:
            if args.bestFit: r = (0.99, 1.01)
            else: r = (0., 2.)

            profiledLoglikelihoodFit = ProfiledLoglikelihoodFit( cardFilePath )
            profiledLoglikelihoodFit.make_workspace(rmin=r[0], rmax=r[1])
            nll = profiledLoglikelihoodFit.likelihoodTest()
            profiledLoglikelihoodFit.cleanup(removeFiles=args.removeCardFiles)
            del profiledLoglikelihoodFit

        logger.info( "NLL: %f", nll)
        ROOT.gDirectory.Clear()

        # in very large WC regions, the fit fails, not relevant for the interesting regions
        if nll is None or abs(nll) > 10000 or abs(nll) < 1: nll = 999

        del c

        return var, nll

    results = []

    SM = calculation( 0 )

    pool = Pool( processes = args.cores )
    results += pool.map( calculation, xRange )
    pool.close()
    del pool

    with open('tmp/'+filename, 'w') as f:
        for item in [SM]+results:
            f.write( "%s\n" % ','.join( map( str, list(item) ) ) )

else:
    with open('dat/'+filename, 'r') as f:
        data = f.readlines()

    results = []
    for i, line in enumerate(data):
        vals = map( float, line.split('\n')[0].split(',') )
        if i == 0:
            if vals[0] != 0:
                raise ValueError('SM Point in data file is not valid!')
            SM = tuple( vals )
        else: results.append( tuple( vals ) )


#Plot

results = [ (x, 2*(res-SM[1])) for x, res in results]
results.sort( key = lambda res: res[0] )

def toGraph( name, title, data ):
    result  = ROOT.TGraph( len(data) )
    result.SetName( name )
    result.SetTitle( title )
    for j, datapoint in enumerate(data):
        x, val = datapoint
        result.SetPoint(j, x, val)
    c = ROOT.TCanvas()
    result.Draw()
    del c
    #res = ROOT.TGraphDelaunay(result)
    return result

# Plot ranges
ranges = {'cpt':[-1.1,1.1], 'cpQM':[-1.1,1.1], 'ctZ':[-0.6,0.6], 'ctZI':[-0.6,0.6]}

polString = "[0]*x**2+[1]*x**3+[2]*x**4+[3]*x**5+[4]*x**6"
# get TGraph from results data list
xhist = toGraph( args.process, args.process, results )
func  = ROOT.TF1("func",polString,ranges[args.variable][0], ranges[args.variable][1] ) 
xhist.Fit(func,"NO")
x68min = func.GetX( 0.989, ranges[args.variable][0], 0 )
x68max = func.GetX( 0.989, 0, ranges[args.variable][1] )
x95min = func.GetX( 3.84, ranges[args.variable][0], 0 )
x95max = func.GetX( 3.84, 0, ranges[args.variable][1] )

xhist.SetLineWidth(0)

func.SetFillColor(ROOT.kWhite)
func.SetFillStyle(1001)
func.SetLineWidth(3)
func.SetLineColor(ROOT.kBlack)
func.SetNpx(1000)

print args.variable, '68', x68min, x68max
print args.variable, '95', x95min, x95max

ROOT.gStyle.SetPadLeftMargin(0.14)
ROOT.gStyle.SetPadRightMargin(0.1)
ROOT.gStyle.SetPadTopMargin(0.11)

# Plot
cans = ROOT.TCanvas("cans","cans",500,500)

if not None in args.yRange:
    xhist.GetYaxis().SetRangeUser( args.yRange[0], args.yRange[1] )
    xhist.GetXaxis().SetRangeUser( ranges[args.variable][0], ranges[args.variable][1] )

func95 = ROOT.TF1("func95",polString, x95min,x95max ) 
xhist.Fit(func95,"NO")
func95.SetFillColor(ROOT.kOrange+7)
func95.SetFillStyle(1001)
func95.SetLineWidth(0)
func95.SetNpx(1000)

func68 = ROOT.TF1("func68",polString, x68min,x68max ) 
xhist.Fit(func68,"NO")
func68.SetFillColor(ROOT.kSpring-1)
func68.SetFillStyle(1001)
func68.SetLineWidth(0)
func68.SetNpx(1000)

if not None in args.yRange:
    func.GetYaxis().SetRangeUser( args.yRange[0], args.yRange[1] )
    func.GetXaxis().SetRangeUser( ranges[args.variable][0], ranges[args.variable][1] )
    func68.GetYaxis().SetRangeUser( args.yRange[0], args.yRange[1] )
    func68.GetXaxis().SetRangeUser( ranges[args.variable][0], ranges[args.variable][1] )
    func95.GetYaxis().SetRangeUser( args.yRange[0], args.yRange[1] )
    func95.GetXaxis().SetRangeUser( ranges[args.variable][0], ranges[args.variable][1] )

xhist.Draw("ALO")
func.Draw("COSAME")
func95.Draw("FOSAME")
func68.Draw("FOSAME")
#xhist.Draw("LOSAME")
func.Draw("COSAME")

# Redraw axis, otherwise the filled graphes overlay
cans.RedrawAxis()

# dashed line at 1
line5 = ROOT.TLine(ranges[args.variable][0], 0.989, ranges[args.variable][1], 0.989 )
line5.SetLineWidth(1)
line5.SetLineStyle(7)
line5.SetLineColor(ROOT.kBlack)
# dashed line at 4
line6 = ROOT.TLine(ranges[args.variable][0], 3.84, ranges[args.variable][1], 3.84 )
line6.SetLineWidth(1)
line6.SetLineStyle(7)
line6.SetLineColor(ROOT.kBlack)

line5.Draw()
line6.Draw()

cans.Update()
xhist.GetYaxis().SetTitle("-2 #Delta ln L")

leg = ROOT.TLegend(0.3,0.7,0.7,0.87)
leg.SetBorderSize(0)
leg.SetTextSize(0.038)
leg.AddEntry(func,"log-likelihood ratio","l")
leg.AddEntry(func68,"68% CL","f")
leg.AddEntry(func95,"95% CL","f")
leg.Draw()

xTitle = args.variable.replace('c','C_{').replace('p','#phi').replace('M','') + '}' 
if 'I' in xTitle: xTitle = xTitle.replace('I','') + '^{[Im]}'
xhist.GetXaxis().SetTitle( xTitle + ' (#Lambda/TeV)^{2}' )

xhist.GetXaxis().SetTitleFont(42)
xhist.GetYaxis().SetTitleFont(42)
xhist.GetXaxis().SetLabelFont(42)
xhist.GetYaxis().SetLabelFont(42)

xhist.GetXaxis().SetTitleOffset(1.3)
xhist.GetYaxis().SetTitleOffset(1.3)

xhist.GetXaxis().SetTitleSize(0.045)
xhist.GetYaxis().SetTitleSize(0.045)
xhist.GetXaxis().SetLabelSize(0.04)
xhist.GetYaxis().SetLabelSize(0.04)

latex1 = ROOT.TLatex()
latex1.SetNDC()
latex1.SetTextSize(0.04)
latex1.SetTextFont(42)
latex1.SetTextAlign(11)

latex1.DrawLatex(0.03, 0.91, '#bf{CMS Phase-2} #it{Simulation Preliminary}'),
latex1.DrawLatex(0.66, 0.91, '%i ab{}^{-1} (%s TeV)' % (int(args.luminosity/1000.), "14" if args.scale14TeV else "13"))

#latex2 = ROOT.TLatex()
#latex2.SetNDC()
#latex2.SetTextSize(0.04)
#latex2.SetTextFont(42)
#latex2.SetTextAlign(11)

#latex2.DrawLatex(0.15, 0.9, 'with Stat. uncert. only' if args.statOnly else 'with YR18 syst. uncert.' if not args.noExpUnc else 'with Stat. and Theory uncert. only'),

plot_directory_ = os.path.join(\
    plot_directory,
    '%s_%s'%(args.level, args.version),
    args.detector,
    args.sample,
    'backgrounds',
    '1Dnll_small' if args.small else '1Dnll',
    args.selection)

if not os.path.isdir( plot_directory_ ):
    os.makedirs( plot_directory_ )

for e in [".png",".pdf",".root"]:
    cans.Print( plot_directory_ + '/' + '_'.join(['1D',args.variable,'lumi'+str(args.luminosity), "14TeV" if args.scale14TeV else "13TeV", "CMScombine" if args.useCombine else "privateFit", "bestFit" if args.bestFit else "r1", 'statOnly' if args.statOnly else 'fullUnc' if not args.noExpUnc else 'noExpUnc']) + e)
