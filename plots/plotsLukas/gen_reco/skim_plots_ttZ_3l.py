#!/usr/bin/env python
''' Analysis script for standard plots
'''

# Standard imports and batch mode
import ROOT, os, itertools
ROOT.gROOT.SetBatch(True)
from math                                import sqrt, cos, sin, pi, isnan, sinh, cosh
import copy
import imp

# RootTools
from RootTools.core.standard             import *

# TTXPheno
from TTXPheno.Tools.user                 import plot_directory
from TTXPheno.Tools.helpers              import deltaPhi, getCollection, deltaR, nanJet, nanLepton, getObjDict
from TTXPheno.Tools.WeightInfo           import WeightInfo

# Import samples
from TTXPheno.samples.benchmarks         import *

# Import helpers
from plot_helpers                        import *

# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',           action='store',      default='INFO', nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--version',            action='store',      default='v7',   help='Appendix to plot directory')
argParser.add_argument('--sample',             action='store',      default='fwlite_ttZ_ll_LO_order2_15weights_ref', help='Sample name specified in sample/python/benchmarks.py, e.g. fwlite_ttZ_ll_LO_order2_15weights_ref')
argParser.add_argument('--order',              action='store',      default=2,      help='Polynomial order of weight string (e.g. 2)')
argParser.add_argument('--selection',          action='store',      default='lepSel3-onZ-njet3p-nbjet1p-Zpt0', help="Specify cut.")
argParser.add_argument('--small',              action='store_true',                 help='Run only on a small subset of the data?')
argParser.add_argument('--level',              action='store',      default='gen',  nargs='?', choices=['reco', 'gen', 'genLep'], help='Which level of reconstruction? reco, gen, genLep')
argParser.add_argument('--scaleLumi',          action='store_true',                 help='Scale lumi only?')
argParser.add_argument('--reweightPtZToSM',    action='store_true',                 help='Reweight Pt(Z) to the SM for all the signals?')
argParser.add_argument('--parameters',         action='store',      default = ['ctW', '3', 'ctWI', '3', 'ctZ', '3', 'ctZI', '3'], type=str, nargs='+', help = "argument parameters")
argParser.add_argument('--luminosity',         action='store',      default=150,    help='Luminosity for weighting the plots')

args = argParser.parse_args()

# Import additional functions/classes specified for the level of reconstruction
if args.level == 'reco':
    from TTXPheno.Tools.cutInterpreterReco   import cutInterpreter
    from TTXPheno.Tools.objectSelection      import isGoodGenJet       as isGoodJet
    from TTXPheno.Tools.objectSelection      import isGoodGenLepton    as isGoodLepton
    preTag = 'reco'
    tag = 'reco'

else:
    from TTXPheno.Tools.cutInterpreter       import cutInterpreter
    from TTXPheno.Tools.objectSelection      import isGoodRecoJet      as isGoodJet
    from TTXPheno.Tools.objectSelection      import isGoodRecoLepton   as isGoodLepton
    preTag = 'gen'
    tag = 'genLep'

# Logger
import TTXPheno.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(    args.logLevel, logFile = None )
logger_rt = logger_rt.get_logger( args.logLevel, logFile = None )

# Make subdirectory
subDirectory = []
if args.scaleLumi: subDirectory.append("shape")
else:              subDirectory.append("lumi")
if args.reweightPtZToSM: subDirectory.append("reweightPtZToSM")
if args.small:     subDirectory.append("small")
subDirectory = '_'.join( subDirectory )

# Import samples
sample_file = "$CMSSW_BASE/python/TTXPheno/samples/benchmarks.py"
samples = imp.load_source( "samples", os.path.expandvars( sample_file ) )
sample = getattr( samples, args.sample )

# Scale the plots with number of events used (implemented in ref_lumiweight1fb)
event_factor = 1.
if args.small:
    sample.reduceFiles( to = 1 )
    event_factor = sample.nEvents / float(sample.chain.GetEntries())

# Polynomial parametrization
w = WeightInfo(sample.reweight_pkl)
w.set_order(int(args.order))

colors = [ ROOT.kMagenta+1, ROOT.kOrange, ROOT.kBlue, ROOT.kCyan+1, ROOT.kGreen+1, ROOT.kRed, ROOT.kViolet, ROOT.kYellow+2 ]

coeffs = args.parameters[::2]
str_vals = args.parameters[1::2]
vals   = list( map( float, str_vals ) )
params = []
for i_param, (coeff, val, str_val) in enumerate(zip(coeffs, vals, str_vals)):
    params.append( { 
        'legendText': ' '.join([coeff,str_val]),
        'WC'        : { coeff:val },
        'color'     : colors[i_param], 
        })

params.append( {'legendText':'SM', 'WC':{}, 'color':ROOT.kBlack} )

# Make stack 
stack  = Stack(*[ [ sample ] for param in params ] )

# reweighting of pTZ 
if args.reweightPtZToSM:
    varZ = "%sZ_pt"%args.level
    for param in params[::-1]:
        param['ptZ_histo'] = sample.get1DHistoFromDraw(varZ, [20,0,500], selectionString = cutInterpreter.cutString(args.selection), weightString = w.get_weight_string(**param['WC']))
        if param['ptZ_histo'].Integral()>0: param['ptZ_histo'].Scale(1./param['ptZ_histo'].Integral())
        param['ptZ_reweight_histo'] = params[-1]['ptZ_histo'].Clone()
        param['ptZ_reweight_histo'].Divide(param['ptZ_histo'])
        logger.info( 'Made reweighting histogram for ptZ and param-point %r with integral %f', param, param['ptZ_reweight_histo'].Integral())

    def get_reweight( param ):
        histo = param['ptZ_reweight_histo']
        bsm_rw = w.get_weight_func( **param['WC'] )
        def reweight(event, sample):
            i_bin = histo.FindBin(getattr( event, varZ ) )
            return histo.GetBinContent(i_bin)*bsm_rw( event, sample ) * event.ref_lumiweight1fb * args.luminosity * event_factor
        return reweight

    weight = [ [ get_reweight( param ) ] for param in params ]
else:
    def get_reweight( param ):

        def reweight(event, sample):
            return w.get_weight_func( **param['WC'] )( event, sample ) * event.ref_lumiweight1fb * args.luminosity * event_factor

        return reweight

    weight = [ [ get_reweight( param ) ] for param in params ]

def drawObjects( hasData = False ):
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextSize(0.04)
    tex.SetTextAlign(11) # align right
    lines = [
      (0.15, 0.95, 'data' if hasData else "Simulation (%s)"%args.level),
      (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV) Scale %3.2f'% ( args.luminosity, dataMCScale ) ) if hasData else (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV)' % args.luminosity)
    ]
    return [tex.DrawLatex(*l) for l in lines] 

def drawPlots(plots):
  for plot in plots:
    for i_h, h in enumerate(plot.histos):
      h[0].style = styles.lineStyle(params[i_h]['color'])

  for log in [False, True]:
    # Directory structure
    WC_directory = '_'.join(args.parameters).rstrip('0').replace('-','m').replace('.','p') if len(args.parameters)>1 else 'SM'
    plot_directory_ = os.path.join(\
        plot_directory,
        '%s_%s'%(args.level, args.version),
        sample.name, 
        subDirectory, 
        args.selection if args.selection is not None else 'no_selection', 
        WC_directory,
        "log" if log else "lin")

    # plot the legend
    l_plot = copy.deepcopy(plots[0])
    for i_h, h in enumerate(l_plot.histos):
      h[0].legendText = params[i_h]['legendText']
      h[0].style = styles.lineStyle(params[i_h]['color'])
      h[0].Scale(0.)
      h[0].GetXaxis().SetTickLength(0.)
      h[0].GetYaxis().SetTickLength(0.)
      h[0].GetXaxis().SetLabelOffset(999.)
      h[0].GetYaxis().SetLabelOffset(999.)
    l_plot.name = "legend"
    l_plot.texX = ''
    l_plot.texY = ''
    plotting.draw(l_plot,
        plot_directory = plot_directory_,
        ratio = None, #{'yRange':(0.1,1.9)} if not args.noData else None,
        logX = False, logY = log, sorting = True,
        #yRange = (0.03, "auto") if log else (0., "auto"),
        #scaling = {i:0 for i in range(1, len(params))} if args.scaleLumi else {}, #Scale BSM shapes to SM (first in list)
        legend =  ( (0.17,0.9-0.05*sum(map(len, l_plot.histos))/3,1.,0.9), 3),
        drawObjects = drawObjects( ),
        copyIndexPHP = True,
    )

    # plot the plots
    for plot in plots:
      for i_h, h in enumerate(plot.histos):
        h[0].legendText = params[i_h]['legendText']
      if not max(l[0].GetMaximum() for l in plot.histos): continue # Empty plot

      plotting.draw(plot,
	    plot_directory = plot_directory_,
	    ratio = None, #{'yRange':(0.1,1.9)} if not args.noData else None,
	    logX = False, logY = log, sorting = True,
	    yRange = (0.03, "auto") if log else (0., "auto"),
            scaling = {i:(len(params)-1) for i in range(len(params)-1)} if args.scaleLumi else {}, #Scale BSM shapes to SM (last in list)
#	    scaling = {i:0 for i in range(1, len(params))} if args.scaleLumi else {}, #Scale BSM shapes to SM (first in list)
	    legend = ( (0.17,0.9-0.05*sum(map(len, plot.histos))/3,1.,0.9), 3),
	    drawObjects = drawObjects( ),
        copyIndexPHP = True,
      )

# Read variables and sequences
# List of variables where gen is replaced by reco for reco
read_variables_gen = [
    "ref_lumiweight1fb/F",
    "genMet_pt/F", "genMet_phi/F",

    "ngenJet/I", "genJet[pt/F,eta/F,phi/F]",
    "ngenLep/I", "genLep[pt/F,eta/F,phi/F,pdgId/I]",

    "genZ_pt/F", "genZ_eta/F", "genZ_phi/F", "genZ_mass/F", "genZ_cosThetaStar/F",

    "genBj0_pt/F", "genBj0_phi/F", "genBj0_eta/F",
    "genBj1_pt/F", "genBj1_phi/F", "genBj1_eta/F",

    "genBjLeadlep_index/I", "genBjLeadhad_index/I",
    "genBjNonZlep_index/I", "genBjNonZhad_index/I",
]

# List of variables where genLep is replaced by reco for reco
read_variables_genLep = [
    "genLepZ_pt/F", "genLepZ_eta/F", "genLepZ_phi/F", "genLepZ_mass/F", "genLepZ_cosThetaStar/F",
    "genLepZ_lldPhi/F", "genLepZ_lldR/F",

    "genLepZ_l1_index/I", "genLepZ_l2_index/I",
    "genLepNonZ_l1_index/I", "genLepNonZ_l2_index/I",
]

if args.level == 'reco':
    read_variables = [ variable.replace('gen', 'reco') for variable in read_variables_gen ]
    read_variables += [ variable.replace('genLep', 'reco') for variable in read_variables_genLep ]
else:
    read_variables = read_variables_gen + read_variables_genLep

read_variables.append( VectorTreeVariable.fromString('p[C/F]', nMax=2000) )

logger.info( "Translating cut %s to %s", args.selection, cutInterpreter.cutString(args.selection) )
sample.setSelectionString( cutInterpreter.cutString(args.selection) )
sample.style = styles.lineStyle(ROOT.kBlue)

#sequence functions
sequence = []

def makeJets( event, sample ):
    ''' Add a list of filtered jets to the event (full list is required for lepton cross cleaning)
    '''
    # Define leptonic b-jets
    event.bj0lep = getObjDict( event, '%sJet_'%preTag, ['pt', 'eta', 'phi'], getattr( event, '%sBjLeadlep_index'%preTag ) ) if getattr( event, '%sBjLeadlep_index'%preTag ) >= 0 else nanLepton()
    event.bj0had = getObjDict( event, '%sJet_'%preTag, ['pt', 'eta', 'phi'], getattr( event, '%sBjLeadhad_index'%preTag ) ) if getattr( event, '%sBjLeadhad_index'%preTag ) >= 0 else nanLepton()

    # Define non-Z b-jets
    event.bjNonZlep = getObjDict( event, '%sJet_'%preTag, ['pt', 'eta', 'phi'], getattr( event, '%sBjNonZlep_index'%preTag ) ) if getattr( event, '%sBjNonZlep_index'%preTag ) >= 0 else nanLepton()
    event.bjNonZhad = getObjDict( event, '%sJet_'%preTag, ['pt', 'eta', 'phi'], getattr( event, '%sBjNonZhad_index'%preTag ) ) if getattr( event, '%sBjNonZhad_index'%preTag ) >= 0 else nanLepton()

    # get (second) hardest bjets
    event.bj0 = {'pt':getattr( event, '%sBj0_pt'%preTag ), 'phi':getattr( event, '%sBj0_phi'%preTag ), 'eta':getattr( event, '%sBj0_eta'%preTag )}
    event.bj1 = {'pt':getattr( event, '%sBj1_pt'%preTag ), 'phi':getattr( event, '%sBj1_phi'%preTag ), 'eta':getattr( event, '%sBj1_eta'%preTag )}

    # Add extra vectors
    for p in [event.bj0, event.bj1, event.bj0lep, event.bj0had, event.bjNonZlep, event.bjNonZhad]:
        addTransverseVector( p )
        addTLorentzVector( p )

    # selection checks
    event.foundBjNonZlep = getattr( event, '%sBjNonZlep_index'%preTag ) >= 0 and isGoodJet( event.bjNonZlep )
#    event.foundBjNonZhad = getattr( event, '%sBjNonZhad_index'%preTag ) >= 0 and isGoodJet( event.bjNonZhad )
#    event.foundBj0lep    = getattr( event, '%sBjLeadlep_index'%preTag ) >= 0 and isGoodJet( event.bj0lep )
#    event.foundBj0had    = getattr( event, '%sBjLeadhad_index'%preTag ) >= 0 and isGoodJet( event.bj0had )

    # choose your selection on b-jets
    event.passing_bjets = event.foundBjNonZlep

sequence.append( makeJets )

def makeMET( event, sample ):
    ''' Make a MET vector to facilitate further calculations
    '''
    event.MET = {'pt':getattr(event, '%sMet_pt'%preTag), 'phi':getattr(event, '%sMet_phi'%preTag)}
    addTransverseVector( event.MET )

sequence.append( makeMET )

def makeZ( event, sample ):
    ''' Make a Z vector to facilitate further calculations (either recoZ, genLepZ or genZ)
    '''
    event.Z_unitVec2D = UnitVectorT2( getattr( event, '%sZ_phi'%args.level ) )
    event.Z_vec4D     = ROOT.TLorentzVector()
    event.Z_vec4D.SetPtEtaPhiM( getattr( event, '%sZ_pt'%args.level ), getattr( event, '%sZ_eta'%args.level ), getattr( event, '%sZ_phi'%args.level ), getattr( event, '%sZ_mass'%args.level ) )
    event.Z_unitVec3D = event.Z_vec4D.Vect().Unit()

sequence.append( makeZ )

def makeLeps( event, sample ):
    ''' Add important leptons (no full list of leptons is required for now)
    '''

    # Define Z leptons
    event.Z_l0 = getObjDict( event, '%sLep_'%preTag, ['pt', 'eta', 'phi', 'pdgId'], getattr( event, '%sZ_l1_index'%tag ) ) if getattr( event, '%sZ_l1_index'%tag ) >= 0 else nanLepton()
    event.Z_l1 = getObjDict( event, '%sLep_'%preTag, ['pt', 'eta', 'phi', 'pdgId'], getattr( event, '%sZ_l2_index'%tag ) ) if getattr( event, '%sZ_l2_index'%tag ) >= 0 else nanLepton()

    # Define non-Z leptons
    event.NonZ_l0 = getObjDict( event, '%sLep_'%preTag, ['pt', 'eta', 'phi', 'pdgId'], getattr( event, '%sNonZ_l1_index'%tag ) ) if getattr( event, '%sNonZ_l1_index'%tag ) >= 0 else nanLepton()
    event.NonZ_l1 = getObjDict( event, '%sLep_'%preTag, ['pt', 'eta', 'phi', 'pdgId'], getattr( event, '%sNonZ_l2_index'%tag ) ) if getattr( event, '%sNonZ_l2_index'%tag ) >= 0 else nanLepton()

    # Add extra vectors
    for p in [ event.Z_l0, event.Z_l1, event.NonZ_l0, event.NonZ_l1 ]:
        addTransverseVector( p )
        addTLorentzVector( p )

    # We may loose some events by cross-cleaning or by thresholds.
    event.foundZl0     = getattr( event, '%sZ_l1_index'%tag ) >= 0 and isGoodLepton( event.Z_l0 )
    event.foundZl1     = getattr( event, '%sZ_l2_index'%tag ) >= 0 and isGoodLepton( event.Z_l1 )
    event.foundZ       = event.Z_l0['pdgId'] * event.Z_l1['pdgId'] < 0 and abs(event.Z_l0['pdgId']) == abs(event.Z_l1['pdgId'])
    event.found3lep    = getattr( event, '%sNonZ_l1_index'%tag ) >= 0 and isGoodLepton( event.NonZ_l0 )

    # choose your selection on leptons
    event.passing_leptons = event.found3lep and event.foundZl0 and event.foundZl1 and event.foundZ

sequence.append( makeLeps )

def makeObservables( event, sample):
    ''' Compute all relevant observables
    '''
    # double b kinematic
    event.bbdPhi = deltaPhi( event.bj0['phi'], event.bj1['phi'] )
    event.bbdR   = deltaR( event.bj0, event.bj1 )

    # Make leptonic W, Nan if any component is Nan
    event.Wlep_vec2D = event.MET['vec2D'] + event.NonZ_l0['vec2D']
    # Lp, Nan if any component is Nan
    event.Wlep_Lp = ( event.Wlep_vec2D*event.NonZ_l0['vec2D'] ) / ( event.Wlep_vec2D*event.Wlep_vec2D )
    # classic MT, Nan if any component is Nan
    event.Wlep_MT = sqrt( MTSquared( event.MET, event.NonZ_l0 ) )

    # blep+lep subsystem, Nan if len(event.lepsNotFromZ == 0)
    event.bleplep_vec2D = event.bjNonZlep['vec2D'] + event.NonZ_l0['vec2D']
    event.bleplep_vec4D = event.bjNonZlep['vec4D'] + event.NonZ_l0['vec4D']

    # transverse mass of top, Nan if len(event.lepsNotFromZ == 0)
    event.t_MT = sqrt( MTSquared( event.NonZ_l0, event.MET ) + MSquared( event.NonZ_l0, event.bjNonZlep ) + MTSquared( event.MET, event.bjNonZlep ) )
    event.t_vec2D = event.NonZ_l0['vec2D'] + event.MET['vec2D'] + event.bjNonZlep['vec2D']

    # signed lepton pt, Nan if len(event.lepsNotFromZ == 0)
    event.getnonZlepchargept = event.NonZ_l0['pt'] if event.NonZ_l0['pdgId']>0 else -event.NonZ_l0['pt']

    # choose your final selection
    event.passing_checks = event.passing_leptons and event.passing_bjets

sequence.append( makeObservables )

# Use some defaults
Plot.setDefaults(stack = stack, weight = weight, addOverFlowBin=None)
  
if args.scaleLumi: y_label = 'norm. diff. xsec'
else: y_label = 'diff. x-sec'

plots = []

plots.append(Plot( name = "Z_pt",
  texX = 'p_{T}(Z) [GeV]', texY = y_label,
  attribute = lambda event, sample: getattr( event, '%sZ_pt'%args.level ) if event.passing_checks else float('nan'),
  binning=[20,0,500],
))

plots.append(Plot( name = "Z_mass",
  texX = 'm(ll) [GeV]', texY = y_label,
  attribute = lambda event, sample: getattr( event, '%sZ_mass'%args.level ) if event.passing_checks else float('nan'),
  binning=[20,70,110],
))

plots.append(Plot( name = 'Z_phi',
  texX = '#phi(Z) [GeV]', texY = y_label,
  attribute = lambda event, sample: getattr( event, '%sZ_phi'%args.level ) if event.passing_checks else float('nan'),
  binning=[20,-pi,pi],
))

plots.append(Plot( name = 'Z_eta',
  texX = '#eta(Z) [GeV]', texY = y_label,
  attribute = lambda event, sample: getattr( event, '%sZ_eta'%args.level ) if event.passing_checks else float('nan'),
  binning=[20,-3,3],
))

plots.append(Plot( name = "Z_cosThetaStar",
  texX = 'cos(#theta*)', texY = y_label,
  attribute = lambda event, sample: getattr( event, '%sZ_cosThetaStar'%args.level ) if event.passing_checks else float('nan'),
  binning=[20,-1.2,1.2],
))

plots.append(Plot( name = "b0_pt",
  texX = 'p_{T}(b_{0}) [GeV]', texY = y_label,
  attribute = lambda event, sample: event.bj0['pt'] if event.passing_checks else float('nan'),
  binning=[20,0,400],
))

plots.append(Plot( name = "b1_pt",
  texX = 'p_{T}(b_{1}) [GeV]', texY = y_label,
  attribute = lambda event, sample: event.bj1['pt'] if event.passing_checks else float('nan'),
  binning=[20,0,400],
))

plots.append(Plot( name = "b0_eta",
  texX = '#eta(b_{0})', texY = y_label,
  attribute = lambda event, sample: event.bj0['eta'] if event.passing_checks else float('nan'),
  binning=[20,-3,3],
))

plots.append(Plot( name = "b1_eta",
  texX = '#eta(b_{1})', texY = y_label,
  attribute = lambda event, sample: event.bj1['eta'] if event.passing_checks else float('nan'),
  binning=[20,-3,3],
))

plots.append(Plot( name = "b0_phi",
  texX = '#phi(b_{0})', texY = y_label,
  attribute = lambda event, sample: event.bj0['phi'] if event.passing_checks else float('nan'),
  binning=[20,pi,pi],
))

plots.append(Plot( name = "b1_phi",
  texX = '#phi(b_{1})', texY = y_label,
  attribute = lambda event, sample: event.bj1['phi'] if event.passing_checks else float('nan'),
  binning=[20,pi,pi],
))

plots.append(Plot( name = 'deltaPhi_bb',
  texX = '#Delta#phi(bb)', texY = y_label,
  attribute = lambda event, sample: event.bbdPhi if event.passing_checks else float('nan'),
  binning=[20,0,pi],
))

plots.append(Plot( name = 'deltaR_bb',
  texX = '#DeltaR(bb)', texY = y_label,
  attribute = lambda event, sample: event.bbdR if event.passing_checks else float('nan'),
  binning=[20,0,6],
))

# for gen use genLepZ_lldPhi
plots.append(Plot( name = 'Z_deltaPhi_ll',
  texX = '#Delta#phi(ll)', texY = y_label,
  attribute = lambda event, sample: getattr( event, '%sZ_lldPhi'%tag ) if event.passing_checks else float('nan'),
  binning=[20,0,pi],
))
    
# for gen use genLepZ_lldR
plots.append(Plot( name = 'Z_deltaR_ll',
  texX = '#DeltaR(ll)', texY = y_label,
  attribute = lambda event, sample: getattr( event, '%sZ_lldR'%tag ) if event.passing_checks else float('nan'),
  binning=[20,0,4],
))
    
plots.append(Plot( name = 'Met_pt',
  texX = 'E_{T}^{miss} [GeV]', texY = y_label,
  attribute = lambda event, sample: event.MET['pt'] if event.passing_checks else float('nan'),
  binning=[20,0,200],
))

plots.append(Plot( name	= 'Met_phi',
  texX = '#phi(E_{T}^{miss})', texY = y_label,
  attribute = lambda event, sample: event.MET['phi'] if event.passing_checks else float('nan'),
  binning=[20,-pi,pi],
))

plots.append(Plot( name = 'W_pt',
  texX = 'p_{T}(W_{lep}) [GeV]', texY = y_label,
  attribute = lambda event, sample: event.Wlep_vec2D.Mod() if event.passing_checks else float('nan'),
  binning=[20,0,400],
))

plots.append(Plot( name = 'W_Lp',
  texX = 'L_{p} from W_{lep}', texY = y_label,
  attribute = lambda event, sample: event.Wlep_Lp if event.passing_checks else float('nan'),
  binning=[20,-3,3],
))

plots.append(Plot( name = 'bleplep_dot_nZ_2D',
  texX = 'p_{T}(b_{lep} + l) [GeV] #upoint n(Z) (2D)', texY = y_label,
  attribute = lambda event, sample: event.bleplep_vec2D*event.Z_unitVec2D if event.passing_checks else float('nan'),
  binning=[20,-400,400],
))

plots.append(Plot( name = 'bleplep_dot_nZ_3D',
  texX = 'p_{T}(b_{lep} + l) [GeV] #upoint n(Z) (3D)', texY = y_label,
  attribute = lambda event, sample: event.bleplep_vec4D.Vect()*event.Z_unitVec3D if event.passing_checks else float('nan'),
  binning=[20,-400,400],
))

plots.append(Plot( name = 'top_dot_nZ',
  texX = 'p_{T}(t_{lep}) [GeV] #upoint n(Z)', texY = y_label,
  attribute = lambda event, sample: event.t_vec2D*event.Z_unitVec2D if event.passing_checks else float('nan'),
  binning=[20,-400,400],
))

plots.append(Plot( name = 'top_lep_pt',
  texX = 'p_{T}(t_{lep}) [GeV]', texY = y_label,
  attribute = lambda event, sample: event.t_vec2D.Mod() if event.passing_checks else float('nan'),
  binning=[20,0,400],
))

plots.append(Plot( name = 'lnonZ_phi',
  texX = '#phi(l^{non-Z}) [GeV]', texY = y_label,
  attribute = lambda event, sample: event.NonZ_l0['phi'] if event.passing_checks else float('nan'),
  binning=[20,-pi,pi],
))

plots.append(Plot( name = 'lnonZ_eta',
  texX = '#eta(l^{non-Z}) [GeV]', texY = y_label,
  attribute = lambda event, sample: event.NonZ_l0['eta'] if event.passing_checks else float('nan'),
  binning=[20,-3,3],
))

plots.append(Plot( name = 'lnonZ_pt',
  texX = 'p_{T}(l^{non-Z}) [GeV]', texY = y_label,
  attribute = lambda event, sample: event.NonZ_l0['pt'] if event.passing_checks else float('nan'),
  binning=[20,0,200],
))

plots.append(Plot( name = 'lnonZ_pt_charge',
  texX = 'p_{T}(l^{non-Z}) [GeV] signed with lepton charge', texY = y_label,
  attribute = lambda event, sample: event.getnonZlepchargept if event.passing_checks else float('nan'),
  binning=[20,-200,200],
))

plots.append(Plot( name = 'mT_W',
  texX = 'm_{T}(W_{lep}) [GeV]', texY = y_label,
  attribute = lambda event, sample: event.Wlep_MT if event.passing_checks else float('nan'),
  binning=[20,0,150],
))

plots.append(Plot( name = 'mT_t',
  texX = 'm_{T}(t_{lep}) [GeV]', texY = y_label,
  attribute = lambda event, sample: event.t_MT if event.passing_checks else float('nan'),
  binning=[20,0,300],
))

plotting.fill(plots, read_variables = read_variables, sequence = sequence, max_events = -1 if args.small else -1)

drawPlots(plots)

#logger.info( "Done with prefix %s and selectionString %s", args.selection, cutInterpreter.cutString(args.selection) )

