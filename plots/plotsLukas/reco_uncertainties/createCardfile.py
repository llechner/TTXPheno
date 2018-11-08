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
argParser.add_argument('--sample',             action='store',     default='fwlite_ttZ_ll_LO_order2_15weights_ref', help='Sample name specified in sample/python/benchmarks.py, e.g. fwlite_ttZ_ll_LO_order2_15weights_ref')
argParser.add_argument('--order',              action='store',     default=2, help='Polynomial order of weight string (e.g. 2)')
argParser.add_argument('--selection',          action='store',     default='lepSel3-onZ-njet3p-nbjet1p-Zpt0-leptonIso3', help="Specify cut.")
argParser.add_argument('--small',              action='store_true', help='Run only on a small subset of the data?')
argParser.add_argument('--luminosity',         action='store',     default=3000, help='Luminosity for weighting the plots')
argParser.add_argument('--detector',           action='store',     default='phase2_CMS', nargs='?', choices=['CMS', 'ATLAS', 'phase2_CMS'], help='Which Delphes detector simulation?')
argParser.add_argument('--scale14TeV',         action='store_true', help='scale 13 TeV cross-sections to 14 Tev?')
argParser.add_argument('--additionalCardFile', action='store',     default='TopEFTCardFile.txt', help='Cardfile where additional uncertainties are taken from')
#argParser.add_argument('--addUncertainties',   action='store',     default = ['trigger_2016','scale','scale_sig','PDF','PartonShower','WZ_xsec', 'ttX'], type=str, help = "add additional uncertainties from cardFile")
argParser.add_argument('--addUncertainties',   action='store',     default = [], type=str, help = "add additional uncertainties from cardFile")
argParser.add_argument('--addBinNumberShift',  action='store',     default = 0, type=int, help = "which bin number does the region start in the additional card file?")
argParser.add_argument('--uncertaintyScale',   action='store',     default = 0.5, type=float, help = "scale factor for additional uncertainties")

args = argParser.parse_args()

if args.small: from TTXPheno.Analysis.regions import recottZRegionsSmall as regions
#else:          from TTXPheno.Analysis.regions import recottZRegionsCosPTZ200 as regions
else:          from TTXPheno.Analysis.regions import recottZRegionsPTZ as regions

# Import additional functions/classes specified for the level of reconstruction
from TTXPheno.Tools.cutInterpreterReco import cutInterpreter

#do the calculation
if True:
        # Import samples
        sample_file     = "$CMSSW_BASE/python/TTXPheno/samples/benchmarks.py"
        loadedSamples   = imp.load_source( "samples", os.path.expandvars( sample_file ) )

        ttXSample       = getattr( loadedSamples, args.sample + '_%s' %args.detector )
        WZSample        = getattr( loadedSamples, 'fwlite_WZ_lep_LO_order2_15weights_%s' %args.detector )
        tWZSample       = getattr( loadedSamples, 'fwlite_tWZ_LO_order2_15weights_%s' %args.detector )
        tZqSample       = getattr( loadedSamples, 'fwlite_tZq_LO_order2_15weights_%s' %args.detector )
        ttgammaSample   = getattr( loadedSamples, 'fwlite_ttgamma_bg_LO_order2_15weights_%s' %args.detector )

        bg = [ WZSample, tWZSample, tZqSample, ttgammaSample ]

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

        for i_region, region in enumerate(regions):
                # compute signal yield for this region (this is the final code)

                logger.info( "At region %s", region )

                # ttX SM
                ttX_coeffList[region]  = ttXSample.weightInfo.getCoeffListFromDraw( ttXSample, selectionString = region.cutString() )
                ttX_SM_rate[region]    = ttXSample.weightInfo.get_weight_yield( ttX_coeffList[region] )

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

                # Our expected observation :-)
                observation[region] = int( round( sum( background_rate[region].values() ) + ttX_SM_rate[region] ) )


        # Write temporary card file
        from TTXPheno.Tools.cardFileWriter import cardFileWriter

        def createCardfile( ):
        #def calculation( var1, var2 ):

#            cardname = 'Cos_SM_cardfile'
            cardname = 'PTZ_SM_cardfile'
            cardFilePath = os.path.join( cardfileLocation, cardname + '.txt' )
            kwargs = {}
    
            c = cardFileWriter.cardFileWriter()

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
                c.addBin(bin_name, ['_'.join(s.name.split('_')[1:3]) for s in bg], nice_name)

                c.specifyObservation( bin_name, observation[region] )

                c.specifyExpectation( bin_name, 'signal', signal_rate[region] )

                c.specifyFlatUncertainty( 'lumi', 1.01 )
                c.specifyUncertainty( 'JES',        bin_name, 'signal', signal_jes_uncertainty[region]        )
                c.specifyUncertainty( 'btagging',   bin_name, 'signal', signal_btagging_uncertainty[region]   )
                c.specifyUncertainty( 'mistagging', bin_name, 'signal', signal_mistagging_uncertainty[region] )
                c.specifyUncertainty( 'muonId',     bin_name, 'signal', signal_muonId_uncertainty[region]     )
                c.specifyUncertainty( 'electronId', bin_name, 'signal', signal_electronId_uncertainty[region] )

                for unc in args.addUncertainties:
                    c.specifyUncertainty( unc,      bin_name, 'signal', 1+(getUncertaintyValue( args.additionalCardFile, args.addBinNumberShift + i_region, 'signal', unc )-1)*args.uncertaintyScale )

                #c.specifyExpectation( bin_name, 'ttX_SM', ttX_SM_rate[region] )
                #c.specifyUncertainty( 'JES', bin_name, 'ttX_SM', ttX_SM_jes_uncertainty[region])
                #c.specifyUncertainty( 'btagging',bin_name, 'ttX_SM', ttX_SM_btagging_uncertainty[region])

                for background in bg:
                    c.specifyExpectation( bin_name, '_'.join( background.name.split('_')[1:3] ), background_rate[region][background.name] )
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


        createCardfile()

