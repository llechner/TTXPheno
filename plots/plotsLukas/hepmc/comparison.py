#!/usr/bin/env python
''' Analysis script for standard plots
'''

# Standard imports
import ROOT, os, imp, sys, copy
ROOT.gROOT.SetBatch(True)
import itertools
import pickle
from math                               import isnan, ceil, pi

# TTXPheno
from TTXPheno.Tools.helpers              import getCollection, getObjDict

# RootTools
from RootTools.core.standard            import *

# Internal Imports
from TTXPheno.Tools.user                import plot_directory
from TTXPheno.Tools.hepMCCutInterpreter import cutInterpreter

# Default Parameter
loggerChoices = ['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET']

# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',           action='store',      default='INFO', nargs='?', choices=loggerChoices,                                help="Log level for logging")
argParser.add_argument('--selection',          action='store',      default='all')#lepSel3-offZ-nJet2p-nBJet1p')
argParser.add_argument('--version',            action='store',      default='v0')
argParser.add_argument('--sample',             action='store',      default='tt', choices = ["tt", "ttZ"])
argParser.add_argument('--pdf',                action='store',      default='1d0')
argParser.add_argument('--small',              action='store_true',                                                                                  help='Run only on a small subset of the data?', )
argParser.add_argument('--normalize',          action='store_true', default=False,                                                                   help="Normalize yields" )
args = argParser.parse_args()

# Logger
import TTXPheno.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

# Samples
from TTXPheno.samples.hepmc_samples      import *
hepSample = ttbarZ if args.sample == "ttZ" else ttbar
hepSample.root_samples_dict = { name:sample for name, sample in hepSample.root_samples_dict.iteritems() if args.pdf in name or name == "PP"}

sample_directory = hepSample.name
if args.small:     sample_directory += "_small"
if args.normalize: sample_directory += "_normalize"

def drawObjects( lumi_scale ):
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextSize(0.04)
    tex.SetTextAlign(11) # align right
    lines = [
      (0.15, 0.95, 'Higgs-PDF Simulation (%s)'%hepSample.name.replace("tbar","#bar{t}")), 
      (0.65, 0.95, '%3.1f fb{}^{-1} (13 TeV)' %lumi_scale)
    ]
    return [tex.DrawLatex(*l) for l in lines] 

if args.normalize:
    scaling = { i:0 for i, _ in enumerate(hepSample.root_samples_dict.keys()) }

# Plotting
def drawPlots( plots, mode ):
    for log in [False, True]:
        plot_directory_ = os.path.join( plot_directory, 'hepmcPlots_%s'%args.version, sample_directory, args.selection, args.pdf, "log" if log else "lin" )

        for plot in plots:
            if not max(l[0].GetMaximum() for l in plot.histos): 
                continue # Empty plot
            postFix = ""# (legacy)"
            extensions_ = ["pdf", "png", "root"] if mode in ['all', 'SF', 'mue'] else ['png']

            plotting.draw( plot,
	                       plot_directory = plot_directory_,
                           extensions = extensions_,
                           ratio = {'yRange': (0.2, 1.8), 'histos':[(i,0) for i, _ in enumerate(hepSample.root_samples_dict.keys())], 'texY':'Ratio'},
#	                       ratio = None,
	                       logX = False, logY = log, sorting = True,
	                       yRange = (0.03, "auto") if log else (0.001, "auto"),
	                       scaling = scaling if args.normalize else {},
	                       legend = [ (0.25,0.88-0.02*sum(map(len, plot.histos)),0.9,0.88), 2],
	                       drawObjects = drawObjects( lumi_scale ) if not args.normalize else drawObjects( lumi_scale ),
                           copyIndexPHP = True,
                         )

def getYieldPlot( index ):
    return Plot(
                name      = 'yield',
                texX      = 'yield',
                texY      = 'Number of Events',
                attribute = lambda event, sample: 0.5 + index,
                binning   = [ 3, 0, 3 ],
                )

recoJetVarString      = "pt/F,eta/F,phi/F,bTag/F,bTagPhys/I,nCharged/I,nNeutrals/I,pt_JEC_up/F,pt_JEC_up/F,bTag_loose/I,bTag_medium/I,bTag_tight/I,bTag_looswMTD/I,bTag_mediumMTD/I,bTag_tightMTD/I"
recoJetVars           = [ item.split("/")[0] for item in recoJetVarString.split(",") ]

recoTopVarString      = "pt/F"#,eta/F,phi/F,mass/F"
recoTopVars           = [ item.split("/")[0] for item in recoTopVarString.split(",") ]

recoLeptonVarString   = "pt/F,eta/F,phi/F,pdgId/I,isolationVar/F,isolationVarRhoCorr/F,sumPtCharged/F,sumPtNeutral/F,sumPtChargedPU/F,sumPt/F,ehadOverEem/F"
recoLeptonVars        = [ item.split("/")[0] for item in recoLeptonVarString.split(",") ]

recoPhotonVarString   = "pt/F,phi/F,eta/F,isolationVar/F,isolationVarRhoCorr/F,sumPtCharged/F,sumPtNeutral/F,sumPtChargedPU/F,sumPt/F,ehadOverEem/F,genIndex/I,minLeptonDR/F,minLeptonPt/F,minJetDR/F"
recoPhotonVars        = [ item.split("/")[0] for item in recoPhotonVarString.split(",") ]



# Read variables and sequences
read_variables  = ["lumiweight1fb/F",

                   "nBTag/I",
                   "nBTag_JEC_down/I","nBTag_JEC_up/I",
                   "nBTag_loose/I","nBTag_medium/I","nBTag_tight/I","nBTag_looswMTD/I","nBTag_mediumMTD/I","nBTag_tightMTD/I",

                   "recoBjNonZlep_index/I",
                   "recoBjNonZhad_index/I",
                   "recoBjLeadlep_index/I",
                   "recoBjLeadhad_index/I",

                   "nrecoPhoton/I",
                   "recoPhoton[%s]"   %recoPhotonVarString,

                   "nrecoJets_JEC_down/I","nrecoJets_JEC_up/I",
                   "nrecoJet/I",
                   "recoJet[%s]"      %recoJetVarString,

                   "nrecoLep/I",
                   "recoLep[%s]"   %recoLeptonVarString,

                   "recoMet_pt/F", "recoMet_phi/F",

                   "recoZ_l1_index/I",
                   "recoZ_l2_index/I",
                   "recoNonZ_l1_index/I",
                   "recoNonZ_l2_index/I",
                   "recoZ_pt/F",
                   "recoZ_eta/F",
                   "recoZ_phi/F",
                   "recoZ_mass/F",
                   "recoZ_lldPhi/F",
                   "recoZ_lldR/F",
                   "recoZ_cosThetaStar/F",

                   "reweight_BTag_B/F",
                   "reweight_BTag_L/F",
                   "reweight_id_mu/F",
                   "reweight_id_ele/F",

                  ]

read_variables += [ "recoBj0_" + var for var in recoJetVarString.split(",") ]
read_variables += [ "recoBj1_" + var for var in recoJetVarString.split(",") ]


#colors
colors = {'HH':ROOT.kOrange+10, 'HG':ROOT.kViolet+6, 'GH':ROOT.kBlue+2}

def addBTag( event, sample ):
    event.jets = getCollection( event, 'recoJet', [ "bTag" ], 'nrecoJet' )
    event.nBTag = len( filter( lambda j: j["bTag"], event.jets ) )

def makeObservables( event, sample ):
    ''' Compute all relevant observables
    '''
    # double b kinematic
    event.lldPhi = deltaPhi( event.recoLep_phi[0], event.recoLep_phi[1] )
    event.lldR   = deltaR( { "pt":event.recoLep_pt[0], "eta":event.recoLep_eta[0], "phi":event.recoLep_phi[0] }, { "pt":event.recoLep_pt[0], "eta":event.recoLep_eta[0], "phi":event.recoLep_phi[0] } )

def printObjects(event, sample):
    #print "lep pt 0", event.recoLep_pt[0]
    #print "lep pt 1", event.recoLep_pt[1]
    #print "lep pt 2", event.recoLep_pt[3]
    #print "lep eta 0", event.recoLep_eta[0]
    #print "lep eta 1", event.recoLep_eta[1]
    #print "lep eta 2", event.recoLep_eta[2]
    #print "lep pdgid 0", event.recoLep_pdgId[0]
    #print "lep pdgid 1", event.recoLep_pdgId[1]
    #print "lep pdgid 2", event.recoLep_pdgId[2]
    #print "Z mass", event.recoZ_mass
    #print "nLep", event.nrecoLep
    #print "nJet", event.nrecoJet
    print "nBTag med", event.nBTag_medium
    print "nBTag", event.nBTag


# Sequence
sequence = [\
            addBTag,
#            makeObservables,
#            printObjects,
           ]


lumi_scale = 136.6
comparisonSamples = []

# Sample definition
for name, sample in hepSample.root_samples_dict.iteritems():
#    sample.weight         = get_reweight( param, sample )
#    sample.read_variables = read_variables_EFT

    if name == "PP":
        sample.texName = name
        sample.style   = styles.lineStyle( ROOT.kBlack, width=3  )
        comparisonSamples.insert( 0, [sample] )
    else:
        sample.texName = "%s (%s)" %(name.split("_")[1], name.split("_")[0].replace("d","."))
        sample.style   = styles.lineStyle( colors[name.split("_")[1]], width=2, dashed=True  )
        comparisonSamples.append( [sample] )

stack      = Stack( *comparisonSamples )

eventScale = 1.
if args.small:
    for sample in stack.samples:
        sample.normalization=1.
        sample.reduceFiles( factor=10 )
        eventScale = 1./sample.normalization

weight_ = lambda event, sample: event.lumiweight1fb*lumi_scale*eventScale

# Use some defaults (set defaults before you create/import list of Plots!!)
Plot.setDefaults( stack=stack, weight=staticmethod( weight_ ), selectionString=cutInterpreter.cutString( args.selection ) )#, addOverFlowBin='upper' )

def getPlots():
    plotList = []

    plotList.append( Plot( name = "dl_pt",
      texX = 'p_{T}(ll) [GeV]', texY = "Number of Events",
      attribute = lambda event, sample: event.recoZ_pt,
      binning=[20,0,400],
    ) )
        
    plotList.append( Plot( name = 'dl_phi',
      texX = '#phi(ll)', texY = "Number of Events",
      attribute = lambda event, sample: event.recoZ_phi,
      binning=[20,-pi,pi],
    ) )

    plotList.append( Plot( name = 'dl_eta',
      texX = '#eta(ll)', texY = "Number of Events",
      attribute = lambda event, sample: event.recoZ_eta,
      binning=[20,-3,3],
    ) )
    
    plotList.append( Plot( name = "dl_mass",
      texX = 'm(ll) [GeV]', texY = "Number of Events",
      attribute = lambda event, sample: event.recoZ_mass,
      binning=[20,70,110],
    ) )
    
    plotList.append( Plot( name = 'dl_dPhi',
      texX = '#Delta#phi(ll)', texY = "Number of Events",
      attribute = lambda event, sample: event.recoZ_lldPhi,
      binning=[20,0,pi],
    ) )

    plotList.append( Plot( name = 'dl_dR',
      texX = '#DeltaR(ll)', texY = "Number of Events",
      attribute = lambda event, sample: event.recoZ_lldR,
      binning=[20,0,5],
    ) )

    plotList.append( Plot( name = 'dl_cosThetaStar',
      texX = 'cos#theta*(ll)', texY = "Number of Events",
      attribute = lambda event, sample: event.recoZ_cosThetaStar,
      binning=[20,-1,1],
    ) )

    plotList.append( Plot( name = "j0_pt",
      texX = 'p_{T}(lead jet) [GeV]', texY = "Number of Events",
      attribute = lambda event, sample: event.recoJet_pt[0],
      binning=[20,0,600],
    ) )

    plotList.append( Plot( name = "j0_phi",
      texX = '#phi(lead jet)', texY = "Number of Events",
      attribute = lambda event, sample: event.recoJet_phi[0],
      binning=[20,-pi,pi],
    ) )
    
    plotList.append( Plot( name = "j0_eta",
      texX = '#eta(lead jet)', texY = "Number of Events",
      attribute = lambda event, sample: event.recoJet_eta[0],
      binning=[20,-3,3],
    ) )

    plotList.append( Plot( name = "j1_pt",
      texX = 'p_{T}(sub-lead jet) [GeV]', texY = "Number of Events",
      attribute = lambda event, sample: event.recoJet_pt[1],
      binning=[20,0,600],
    ) )

    plotList.append( Plot( name = "j1_phi",
      texX = '#phi(sub-lead jet)', texY = "Number of Events",
      attribute = lambda event, sample: event.recoJet_phi[1],
      binning=[20,-pi,pi],
    ) )
    
    plotList.append( Plot( name = "j1_eta",
      texX = '#eta(sub-lead jet)', texY = "Number of Events",
      attribute = lambda event, sample: event.recoJet_eta[1],
      binning=[20,-3,3],
    ) )

    plotList.append( Plot( name = "b0_pt",
      texX = 'p_{T}(b_{0}) [GeV]', texY = "Number of Events",
      attribute = lambda event, sample: event.recoBj0_pt,
      binning=[20,0,400],
    ) )

    plotList.append( Plot( name = "b0_phi",
      texX = '#phi(b_{0})', texY = "Number of Events",
      attribute = lambda event, sample: event.recoBj0_phi,
      binning=[20,-pi,pi],
    ) )
    
    plotList.append( Plot( name = "b0_eta",
      texX = '#eta(b_{0})', texY = "Number of Events",
      attribute = lambda event, sample: event.recoBj0_eta,
      binning=[20,-3,3],
    ) )

    plotList.append( Plot( name = "b1_pt",
      texX = 'p_{T}(b_{1}) [GeV]', texY = "Number of Events",
      attribute = lambda event, sample: event.recoBj1_pt,
      binning=[20,0,400],
    ) )
    
    plotList.append( Plot( name = "b1_phi",
      texX = '#phi(b_{1})', texY = "Number of Events",
      attribute = lambda event, sample: event.recoBj1_phi,
      binning=[20,-pi,pi],
    ) )
    
    plotList.append( Plot( name = "b1_eta",
      texX = '#eta(b_{1})', texY = "Number of Events",
      attribute = lambda event, sample: event.recoBj1_eta,
      binning=[20,-3,3],
    ) )
    
    plotList.append( Plot( name = 'Met_pt',
      texX = 'E_{T}^{miss} [GeV]', texY = "Number of Events",
      attribute = lambda event, sample: event.recoMet_pt,
      binning=[20,0,400],
    ) )
    
    plotList.append( Plot( name = 'Met_phi',
      texX = '#phi(E_{T}^{miss})', texY = "Number of Events",
      attribute = lambda event, sample: event.recoMet_phi,
      binning=[10,-pi,pi],
    ) )

    plotList.append( Plot( name = 'l0_pt',
      texX = 'p_{T}(lead lep) [GeV]', texY = "Number of Events",
      attribute = lambda event, sample: event.recoLep_pt[0],
      binning=[20,0,300],
    ) )
    
    plotList.append( Plot( name = 'l0_phi',
      texX = '#phi(lead lep)', texY = "Number of Events",
      attribute = lambda event, sample: event.recoLep_phi[0],
      binning=[20,-pi,pi],
    ) )

    plotList.append( Plot( name = 'l0_eta',
      texX = '#eta(lead lep)', texY = "Number of Events",
      attribute = lambda event, sample: event.recoLep_eta[0],
      binning=[20,-3,3],
    ) )
    
    plotList.append( Plot( name = 'l1_pt',
      texX = 'p_{T}(sub-lead lep) [GeV]', texY = "Number of Events",
      attribute = lambda event, sample: event.recoLep_pt[1],
      binning=[20,0,300],
    ) )
    
    plotList.append( Plot( name = 'l1_phi',
      texX = '#phi(sub-lead lep)', texY = "Number of Events",
      attribute = lambda event, sample: event.recoLep_phi[1],
      binning=[20,-pi,pi],
    ) )

    plotList.append( Plot( name = 'l1_eta',
      texX = '#eta(sub-lead lep)', texY = "Number of Events",
      attribute = lambda event, sample: event.recoLep_eta[1],
      binning=[20,-3,3],
    ) )
    
    plotList.append( Plot( name = "nJet",
      texX = 'N_{jets}', texY = "Number of Events",
      attribute = lambda event, sample: event.nrecoJet,
      binning=[10,0,10],
    ) )

    plotList.append( Plot( name = "nBJet",
      texX = 'N_{b-jets}', texY = "Number of Events",
      attribute = lambda event, sample: event.nBTag,
      binning=[5,0,5],
    ) )

    plotList.append( Plot( name = "nLepton",
      texX = 'N_{lep}', texY = "Number of Events",
      attribute = lambda event, sample: event.nrecoLep,
      binning=[4,0,4],
    ) )

    plotList.append( Plot( name = "nPhoton",
      texX = 'N_{#gamma}', texY = "Number of Events",
      attribute = lambda event, sample: event.nrecoPhoton,
      binning=[3,0,3],
    ) )

    return plotList

plots = getPlots()

# Loop over channels
allPlots = {}
#allModes = [ 'all', 'mumumu', 'mumue', 'muee', 'eee' ]
allModes = [ 'all' ]

for index, mode in enumerate( allModes ):
    logger.info( "Computing plots for mode %s", mode )

    # Define 2l selections
    leptonSelection = cutInterpreter.cutString( mode )

    for sample in stack.samples: sample.setSelectionString( [ leptonSelection ] )

    plotting.fill( plots, read_variables=read_variables, sequence=sequence )

    logger.info( "Plotting mode %s", mode )
    allPlots[mode] = copy.deepcopy(plots) # deep copy for creating SF/all plots afterwards!
    drawPlots( allPlots[mode], mode )

