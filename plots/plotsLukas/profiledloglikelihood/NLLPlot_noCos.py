''' Plot script WC parameter LogLikelihood
'''

# Standard imports 
import sys
import ROOT
import imp
import pickle
import ctypes
import numpy as np

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
argParser.add_argument('--variables' ,         action='store',     default = ['ctZ', 'ctZI'], type=str, nargs=2, help = "argument plotting variables")
argParser.add_argument('--binning',            action='store',     default = [1, -2, 2, 1, -2, 2], type=float, nargs=6, help = "argument parameters")
argParser.add_argument('--zRange',             action='store',     default = [None, None], type=float, nargs=2, help = "argument parameters")
argParser.add_argument('--luminosity',         action='store',     default=150, help='Luminosity for weighting the plots')
argParser.add_argument('--scale',              action='store',     default=None, help='Luminosity for weighting the plots')
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
        else:          from TTXPheno.Analysis.regions import recottZRegionsPTZOnly as regions
    elif 'ttgamma' in args.process.split('_'):
        if args.small: from TTXPheno.Analysis.regions import recottgammaRegionsSmall as regions
        else:          from TTXPheno.Analysis.regions import recottgammaRegions as regions

    # Import additional functions/classes specified for the level of reconstruction
    from TTXPheno.Tools.cutInterpreterReco import cutInterpreter

#binning range
binningX = args.binning[:3]
binningY = args.binning[3:]

if binningX[0] > 1:
    xRange = np.linspace( binningX[1], binningX[2], int(binningX[0]), endpoint=False)
    xRange = [ el + 0.5 * ( xRange[1] - xRange[0] ) for el in xRange ]
else:
    xRange = [ 0.5 * ( binningX[1] + binningX[2] ) ]

if binningY[0] > 1:
    yRange = np.linspace( binningY[1], binningY[2], int(binningY[0]), endpoint=False)
    yRange = [ el + 0.5 * ( yRange[1] - yRange[0] ) for el in yRange ]
else:
    yRange = [ 0.5 * ( binningY[1] + binningY[2] ) ]

#save data file
filename = '_'.join( ['nll', args.detector ] + args.sample.split('_')[1:3] + args.variables + map( str, args.binning ) + [ args.selection, str(args.luminosity), "14TeV" if args.scale14TeV else "13TeV" ] ) + '_PTZOnly.data'

#do the calculation
if not os.path.isfile('dat/' + filename) or args.overwrite:

    if not args.fitOnly:

        # Import samples
        sample_file     = "$CMSSW_BASE/python/TTXPheno/samples/benchmarks.py"
        loadedSamples   = imp.load_source( "samples", os.path.expandvars( sample_file ) )

        ttXSample       = getattr( loadedSamples, args.sample + '_%s' %args.detector )
        WZSample        = getattr( loadedSamples, 'fwlite_WZ_lep_LO_order2_15weights_%s' %args.detector )
    #    ttSample        = getattr( loadedSamples, 'fwlite_tt_full_LO_order2_15weights_%s' %args.detector )
    #    tWSample        = getattr( loadedSamples, 'fwlite_tW_LO_order2_15weights_%s' %args.detector )
        tWZSample       = getattr( loadedSamples, 'fwlite_tWZ_LO_order2_15weights_%s' %args.detector )
        tZqSample       = getattr( loadedSamples, 'fwlite_tZq_LO_order2_15weights_%s' %args.detector )
    #    ZgammaSample    = getattr( loadedSamples, 'fwlite_Zgamma_LO_order2_15weights_%s' %args.detector )
        ttgammaSample   = getattr( loadedSamples, 'fwlite_ttgamma_bg_LO_order2_15weights_%s' %args.detector )

        #if args.process.split('_')[0] == 'ttgamma':
        #    ttgammaIsrSample  = copy.deepcopy( ttXSample ) #select ttgamma events with isolated gamma from ISR (cat a2)
        #    ttgammaIsrSample.name = 'fwlite_ttgamma_ISR_LO_order2_15weights_ref'

        if args.process == 'ttZ_3l': bg = [ WZSample, tWZSample, tZqSample, ttgammaSample ]
    #    if args.process == 'ttZ_3l': bg = [ WZSample, tWZSample, tZqSample ]
        elif args.process == 'ttZ_4l': bg = [ WZSample, tWZSample, tZqSample, ttgammaSample ]
        elif args.process == 'ttgamma_1l': bg = [ ttSample, tWSample, tWZSample, tZqSample, ZgammaSample ]
        elif args.process == 'ttgamma_2l': bg = [ ttSample, tWSample, tWZSample, tZqSample, ZgammaSample ]

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

        for var in args.variables:
            if var not in ttXSample.weightInfo.variables and not (args.variables[0] == 'cuB' and args.variables[1] == 'cuW'):
                raise ValueError('Input variable not in gridpack: %s' %var)

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


    def calculation( variables ):
    #def calculation( var1, var2 ):

        if args.variables[0] == 'cuB' and args.variables[1] == 'cuW':
            var1, var2 = variables #cuB cuW
            ctZ, ctW = cuBWtoctWZ( var1, var2 )
            kwargs = { 'ctZ':ctZ, 'ctW':ctW }
        else:
            var1, var2 = variables
            kwargs = { args.variables[0]:var1, args.variables[1]:var2 }

        nameList = args.sample.split('_')[1:3] + args.variables + args.binning + [ args.level, args.version, args.order, args.luminosity, "14TeV" if args.scale14TeV else "13TeV", args.selection, 'small' if args.small else 'full', 'statOnly' if args.statOnly else 'fullUnc' if not args.noExpUnc else 'noExpUnc', var1, var2 ]
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

        return var1, var2, nll

    results = []

    SM = calculation( (0, 0) )

    for varX in xRange:
        # do not run all calc in one pool, memory leak!!!
        pool = Pool( processes = args.cores )
        results += pool.map( calculation, [ (varX, varY) for varY in yRange ] )
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
        if args.scale is not None: vals[2] = vals[2]*float(args.scale)/float(args.luminosity)/2
        if i == 0:
            if vals[0] != 0 or vals[1] != 0:
                raise ValueError('SM Point in data file is not valid!')
            SM = tuple( vals )
        else: results.append( tuple( vals ) )


#Plot

#scale to SM
results.sort( key = lambda res: ( abs(res[0]), abs(res[1]), res[2] ) )
nll_SM = SM[2]

results = [ (x, y, 2*(result - nll_SM)) for x, y, result in results ]

def toGraph2D( name, title, data ):
    result = ROOT.TGraph2D( len(data) )
    debug = ROOT.TGraph()
    result.SetName( name )
    result.SetTitle( title )
    for i, datapoint in enumerate(data):
        x, y, val = datapoint
        result.SetPoint(i, x, y, val)
        debug.SetPoint(i, x, y)
    c = ROOT.TCanvas()
    result.Draw()
    debug.Draw()
    del c
    #res = ROOT.TGraphDelaunay(result)
    return result, debug

#get TGraph2D from results list
a, debug = toGraph2D( args.process, args.process, results )#res_dic)
nxbins   = max(1, min(500, int(binningX[0])*int(args.binMultiplier)))
nybins   = max(1, min(500, int(binningY[0])*int(args.binMultiplier)))

#re-bin
hist = a.GetHistogram().Clone()
a.SetNpx(nxbins)
a.SetNpy(nybins)
hist = a.GetHistogram().Clone()

#smoothing
if args.smooth: hist.Smooth()

cans = ROOT.TCanvas("can_%s"%args.process,"",500,500)

#calculate contour lines (1sigma, 2sigma) for 2D
contours = {'ttZ_3l': [1.515*1.515, 2.486*2.486], 'ttgamma_1l': [1.515*1.515, 2.486*2.486], 'ttgamma_2l': [1.515*1.515, 2.486*2.486]}
if args.contours:
    histsForCont = hist.Clone()
    c_contlist = ((ctypes.c_double)*(len(contours[args.process])))(*contours[args.process])
    histsForCont.SetContour(len(c_contlist),c_contlist)
    histsForCont.Draw("contzlist")
    cans.Update()
    conts = ROOT.gROOT.GetListOfSpecials().FindObject("contours")
    #cont_m2 = conts.At(0).Clone()
    #cont_m1 = conts.At(1).Clone()
    cont_p1 = conts.At(0).Clone()
    cont_p2 = conts.At(1).Clone()

pads = ROOT.TPad("pad_%s"%args.process,"",0.,0.,1.,1.)
pads.SetRightMargin(0.20)
pads.SetLeftMargin(0.14)
pads.SetTopMargin(0.11)
pads.Draw()
pads.cd()

hist.Draw("colz")

#draw contour lines
if args.contours:
    for conts in [cont_p2]:
        for cont in conts:
            cont.SetLineColor(ROOT.kOrange+7)
            cont.SetLineWidth(2)
#            cont.SetLineStyle(7)
            cont.Draw("same")
    for conts in [cont_p1]:
        for cont in conts:
            cont.SetLineColor(ROOT.kSpring-1)
            cont.SetLineWidth(2)
#            cont.SetLineStyle(7)
            cont.Draw("same")


hist.GetZaxis().SetTitle("-2 #Delta ln L")

if not None in args.zRange:
    hist.GetZaxis().SetRangeUser( args.zRange[0], args.zRange[1] )
#    hist.GetXaxis().SetRangeUser( -0.3 , 0.3 )
#    hist.GetYaxis().SetRangeUser( -0.3 , 0.3 )
#    hist.GetXaxis().SetRangeUser( -1 , 1 )
#    hist.GetYaxis().SetRangeUser( -1 , 1 )
#    hist.GetXaxis().SetRangeUser( -8 , 12 )
#    hist.GetYaxis().SetRangeUser( -8 , 12 )


if args.variables[0] == 'cuB' and args.variables[1] == 'cuW':
    hist.GetXaxis().SetTitle('C^{(33)}_{uB} (#Lambda/TeV)^{2}' )
    hist.GetYaxis().SetTitle('C^{(33)}_{uW} (#Lambda/TeV)^{2}' )
else:
    xTitle = args.variables[0].replace('c','C_{').replace('p','#phi').replace('M','') + '}' 
    if 'I' in xTitle: xTitle = xTitle.replace('I','') + '^{[Im]}'
    yTitle = args.variables[1].replace('c','C_{').replace('p','#phi').replace('M','') + '}' 
    if 'I' in yTitle: yTitle = yTitle.replace('I','') + '^{[Im]}'
    hist.GetXaxis().SetTitle( xTitle + ' (#Lambda/TeV)^{2}' )
    hist.GetYaxis().SetTitle( yTitle + ' (#Lambda/TeV)^{2}' )

hist.GetXaxis().SetTitleFont(42)
hist.GetYaxis().SetTitleFont(42)
hist.GetZaxis().SetTitleFont(42)
hist.GetXaxis().SetLabelFont(42)
hist.GetYaxis().SetLabelFont(42)
hist.GetZaxis().SetLabelFont(42)

hist.GetXaxis().SetTitleOffset(1.4)
hist.GetYaxis().SetTitleOffset(1.4)

hist.GetXaxis().SetTitleSize(0.04)
hist.GetYaxis().SetTitleSize(0.04)
hist.GetZaxis().SetTitleSize(0.04)
hist.GetXaxis().SetLabelSize(0.04)
hist.GetYaxis().SetLabelSize(0.04)
hist.GetZaxis().SetLabelSize(0.04)

latex1 = ROOT.TLatex()
latex1.SetNDC()
latex1.SetTextSize(0.04)
latex1.SetTextFont(42)
latex1.SetTextAlign(11)

latex1.DrawLatex(0.15, 0.95, '#bf{CMS} #it{Simulation Preliminary}'),
latex1.DrawLatex(0.55, 0.90, '%i fb{}^{-1} (%s TeV)' % (int(args.luminosity), "14" if args.scale14TeV else "13"))

latex2 = ROOT.TLatex()
latex2.SetNDC()
latex2.SetTextSize(0.04)
latex2.SetTextFont(42)
latex2.SetTextAlign(11)

latex2.DrawLatex(0.15, 0.9, 'with Stat. uncert. only' if args.statOnly else 'with YR18 syst. uncert.' if not args.noExpUnc else 'with Stat. and Theory uncert. only'),

#latex1.DrawLatex(0.15, 0.92, ' '.join(args.process.split('_')[:2]) + ' (' + args.detector + ')')
#latex1.DrawLatex(0.55, 0.92, '%3.1f fb{}^{-1} @ 13 TeV'%(float(args.luminosity) if args.scale is None else float(args.scale)) )

plot_directory_ = os.path.join(\
    plot_directory,
    '%s_%s'%(args.level, args.version),
    args.detector,
    args.sample,
    'backgrounds',
    'nll_small' if args.small else 'nll',
    args.selection)

if not os.path.isdir( plot_directory_ ):
    os.makedirs( plot_directory_ )

for e in [".png",".pdf",".root"]:
    cans.Print( plot_directory_ + '/' + '_'.join(args.variables + ['lumi'+str(args.luminosity) if args.scale is None else 'lumi'+str(args.scale), "14TeV" if args.scale14TeV else "13TeV", "CMScombine" if args.useCombine else "privateFit", "bestFit" if args.bestFit else "r1", 'statOnly' if args.statOnly else 'fullUnc' if not args.noExpUnc else 'noExpUnc']) + e)

