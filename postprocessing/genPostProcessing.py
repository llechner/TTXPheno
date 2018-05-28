#!/usr/bin/env python
''' Make flat ntuple from GEN data tier 
'''
#
# Standard imports and batch mode
#
import ROOT
import os, sys
ROOT.gROOT.SetBatch(True)
import itertools
from math                                import sqrt, cos, sin, pi, acos
import imp

#RootTools
from RootTools.core.standard             import *

#TTXPheno
from TTXPheno.Tools.user                   import skim_output_directory
from TTXPheno.Tools.GenSearch              import GenSearch
from TTXPheno.Tools.helpers                import deltaR2, cosThetaStar
from TTXPheno.Tools.HyperPoly              import HyperPoly
from TTXPheno.Tools.WeightInfo             import WeightInfo
from TTXPheno.Tools.DelphesProducer        import DelphesProducer

#
# Arguments
# 
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',           action='store',      default='INFO',          nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--small',              action='store_true', help='Run only on a small subset of the data?')#, default = True)
argParser.add_argument('--delphes',            action='store_true', help='Run Delphes?')
argParser.add_argument('--overwrite',          action='store_true', help='Overwrite?')#, default = True)
argParser.add_argument('--targetDir',          action='store',      default='v3')
argParser.add_argument('--sample',             action='store',      default='fwlite_ttZ_ll_LO_scan', help="Name of the sample loaded from fwlite_benchmarks. Only if no inputFiles are specified")
argParser.add_argument('--inputFiles',         action='store',      nargs = '*', default=[])
argParser.add_argument('--targetSampleName',   action='store',      default=None, help="Name of the sample in case inputFile are specified. Otherwise ignored")
argParser.add_argument('--nJobs',              action='store',      nargs='?', type=int, default=1,  help="Maximum number of simultaneous jobs.")
argParser.add_argument('--job',                action='store',      nargs='?', type=int, default=0,  help="Run only job i")
argParser.add_argument('--addReweights',       action='store_true',   help="Add reweights?")
argParser.add_argument('--interpolationOrder', action='store',      nargs='?', type=int, default=3,  help="Interpolation order for EFT weights.")
args = argParser.parse_args()

#
# Logger
#
import TTXPheno.Tools.logger as _logger
import RootTools.core.logger as _logger_rt
logger    = _logger.get_logger(   args.logLevel, logFile = None)
logger_rt = _logger_rt.get_logger(args.logLevel, logFile = None)

# Load sample either from 
if len(args.inputFiles)>0:
    logger.info( "Input files found. Ignoring 'sample' argument. Files: %r", args.inputFiles)
    sample = FWLiteSample( args.targetSampleName, args.inputFiles)
else:
    sample_file = "$CMSSW_BASE/python/TTXPheno/samples/fwlite_benchmarks.py"
    samples = imp.load_source( "samples", os.path.expandvars( sample_file ) )
    sample = getattr( samples, args.sample )
    logger.debug( 'Loaded sample %s with %i files.', sample.name, len(sample.files) )

# output directory
output_directory = os.path.join(skim_output_directory, 'gen', args.targetDir, sample.name) 
if not os.path.exists( output_directory ): 
    os.makedirs( output_directory )
    logger.info( "Created output directory %s", output_directory )

maxEvents = -1
if args.small: 
    args.targetDir += "_small"
    maxEvents=10 # Number of files
    sample.files=sample.files[:1]

# Load reweight pickle file if supposed to keep weights. 
extra_variables = []
if args.addReweights:
    weightInfo = WeightInfo( sample.reweight_pkl )
    # Determine coefficients for storing in vector
    # Sort Ids wrt to their position in the card file

    # weights for the ntuple
    rw_vector       = TreeVariable.fromString( "rw[w/F,"+",".join(w+'/F' for w in weightInfo.variables)+"]")
    rw_vector.nMax  = weightInfo.nid
    extra_variables.append(rw_vector)

    # coefficients for the weight parametrization
    param_vector      = TreeVariable.fromString( "p[C/F]" )
    param_vector.nMax = HyperPoly.get_ndof(weightInfo.nvar, args.interpolationOrder)
    hyperPoly         = HyperPoly( args.interpolationOrder )
    extra_variables.append(param_vector)
    extra_variables.append(TreeVariable.fromString( "chi2_ndof/F"))

def interpret_weight(weight_id):
    str_s = weight_id.split('_')
    res={}
    for i in range(len(str_s)/2):
        res[str_s[2*i]] = float(str_s[2*i+1].replace('m','-').replace('p','.'))
    return res

# Run only job number "args.job" from total of "args.nJobs"
if args.nJobs>1:
    n_files_before = len(sample.files)
    sample = sample.split(args.nJobs)[args.job]
    n_files_after  = len(sample.files)
    logger.info( "Running job %i/%i over %i files from a total of %i.", args.job, args.nJobs, n_files_after, n_files_before)

products = {
    'lhe':{'type':'LHEEventProduct', 'label':("externalLHEProducer")},
    'gp':{'type':'vector<reco::GenParticle>', 'label':("genParticles")},
    'genJets':{'type':'vector<reco::GenJet>', 'label':("ak4GenJets")},
    'genMET':{'type':'vector<reco::GenMET>',  'label':("genMetTrue")},
}

def varnames( vec_vars ):
    return [v.split('/')[0] for v in vec_vars.split(',')]

# standard variables
variables  = ["run/I", "lumi/I", "evt/l"]
# MET
variables += ["GenMet_pt/F", "GenMet_phi/F"]

# jet vector
jet_read_vars       =  "pt/F,eta/F,phi/F"
jet_read_varnames   =  varnames( jet_read_vars )
jet_write_vars      = jet_read_vars+',matchBParton/I' 
jet_write_varnames  =  varnames( jet_write_vars )
variables     += ["GenJet[%s]"%jet_write_vars]
# lepton vector 
lep_vars       =  "pt/F,eta/F,phi/F,pdgId/I"
lep_extra_vars =  "motherPdgId/I"
lep_varnames   =  varnames( lep_vars ) 
lep_all_varnames = lep_varnames + varnames(lep_extra_vars)
variables     += ["GenLep[%s]"%(','.join([lep_vars, lep_extra_vars]))]
# top vector
top_vars       =  "pt/F,eta/F,phi/F"
top_varnames   =  varnames( top_vars ) 
variables     += ["top[%s]"%top_vars]

# to be stored for each boson
boson_read_varnames= [ 'pt', 'phi', 'eta', 'mass']
# Z vector
variables     += ["Z_pt/F", "Z_phi/F", "Z_eta/F", "Z_mass/F", "Z_cosThetaStar/F", "Z_daughterPdg/I"]
# W vector
variables     += ["W_pt/F", "W_phi/F", "W_eta/F", "W_mass/F", "W_daughterPdg/I"]
# gamma vector
variables     += ["gamma_pt/F", "gamma_phi/F", "gamma_eta/F", "gamma_mass/F"]
if args.addReweights:
    variables.append('rw_nominal/F')

def fill_vector( event, collection_name, collection_varnames, objects):
    setattr( event, "n"+collection_name, len(objects) )
    for i_obj, obj in enumerate(objects):
        for var in collection_varnames:
            getattr(event, collection_name+"_"+var)[i_obj] = obj[var]

reader = sample.fwliteReader( products = products )

def readReferencePoint():
    # Load ref-point from pkl file
    if len( weightInfo.ref_point.keys() ) > 0:
        # list values of ref-point in right order
        logger.info( "Using reference point from the pkl file!" )
        return [ float( weightInfo.ref_point[var] ) if var in weightInfo.ref_point.keys() else 0 for var in weightInfo.variables ]

    # if no ref point in pkl file, return None
    else:
        logger.info( "No reference point found in pkl file! Continuing without a reference point!")
        return None

def filler( event ):

    event.run, event.lumi, event.evt = reader.evt

    if reader.position % 100==0: logger.info("At event %i/%i", reader.position, reader.nEvents)

    if args.addReweights:
        event.nrw = weightInfo.nid
        lhe_weights = reader.products['lhe'].weights()
        weights      = []
        param_points = []
        for weight in lhe_weights:
            # Store nominal weight (First position!) 
            if weight.id=='rwgt_1': event.rw_nominal = weight.wgt
            if not weight.id in weightInfo.id: continue
            pos = weightInfo.data[weight.id]
            event.rw_w[pos] = weight.wgt
            weights.append( weight.wgt )
            interpreted_weight = interpret_weight(weight.id) 
            for var in weightInfo.variables:
                getattr( event, "rw_"+var )[pos] = interpreted_weight[var]
            # weight data for interpolation
            if not hyperPoly.initialized: param_points.append( tuple(interpreted_weight[var] for var in weightInfo.variables) )

        # Initialize with Reference Point
        if not hyperPoly.initialized: hyperPoly.initialize( param_points, readReferencePoint() )
        coeff = hyperPoly.get_parametrization( weights )

        # = HyperPoly(weight_data, args.interpolationOrder)
        event.np = hyperPoly.ndof
        event.chi2_ndof = hyperPoly.chi2_ndof(coeff, weights)
        #logger.debug( "chi2_ndof %f coeff %r", event.chi2_ndof, coeff )
        logger.debug( "chi2_ndof %f", event.chi2_ndof )
        for n in xrange(hyperPoly.ndof):
            event.p_C[n] = coeff[n]

    # All gen particles
    gp      = reader.products['gp']

    # for searching
    search  = GenSearch( gp )

    # find heavy objects before they decay
    tops = map( lambda t:{var: getattr(t, var)() for var in top_varnames}, filter( lambda p:abs(p.pdgId())==6 and search.isLast(p),  gp) )

    tops.sort( key = lambda p:-p['pt'] )
    fill_vector( event, "top", top_varnames, tops ) 

    # generated Z's
    gen_Zs = filter( lambda p:abs(p.pdgId())==23 and search.isLast(p), gp)
    gen_Zs.sort( key = lambda p: -p.pt() )
    if len(gen_Zs)>0: 
        gen_Z = gen_Zs[0]
        for var in boson_read_varnames:
           setattr( event, "Z_"+var,  getattr(gen_Z, var)() )
    else:
        gen_Z = None
    
    if gen_Z is not None:

        d1, d2 = gen_Z.daughter(0), gen_Z.daughter(1)
        if d1.pdgId()>0: 
            lm, lp = d1, d2
        else:
            lm, lp = d2, d1
        event.Z_daughterPdg = lm.pdgId()
        event.Z_cosThetaStar = cosThetaStar(gen_Z.mass(), gen_Z.pt(), gen_Z.eta(), gen_Z.phi(), lm.pt(), lm.eta(), lm.phi())

    # generated W's
    gen_Ws = filter( lambda p:abs(p.pdgId())==24 and search.isLast(p), gp)
    gen_Ws.sort( key = lambda p: -p.pt() )
    # W can't have a top-mother - We're looking for the extra boson (there is always an extra boson)
    gen_Ws = filter( lambda p: abs(search.ascend(p).mother(0).pdgId())!=6, gen_Ws )
    if len(gen_Ws)>0: 
        gen_W = gen_Ws[0]
        for var in boson_read_varnames:
           setattr( event, "W_"+var,  getattr(gen_W, var)() )
    else:
        gen_W = None
    
    if gen_W is not None:

        d1, d2 = gen_W.daughter(0), gen_W.daughter(1)
        if abs(d1.pdgId()) in [11, 13, 15]: 
            lep, neu = d1, d2
        else:
            lep, neu = d2, d1

        event.W_daughterPdg = lep.pdgId()

    gen_Gammas = filter( lambda p:abs(p.pdgId())==22 and search.isLast(p), gp)
    gen_Gammas.sort( key = lambda p: -p.pt() )
    if len(gen_Gammas)>0: 
        gen_Gamma = gen_Gammas[0]
        for var in boson_read_varnames:
           setattr( event, "gamma_"+var,  getattr(gen_Gamma, var)() )
    else:
        gen_Gamma = None
    
    # find all leptons 
    leptons = [ (search.ascend(l), l) for l in filter( lambda p:abs(p.pdgId()) in [11, 13] and search.isLast(p) and p.pt()>=0,  gp) ]
    leps    = []
    for first, last in leptons:
        mother_pdgId = first.mother(0).pdgId() if first.numberOfMothers()>0 else -1
        leps.append( {var: getattr(last, var)() for var in lep_varnames} )
        leps[-1]['motherPdgId'] = mother_pdgId

    leps.sort( key = lambda p:-p['pt'] )
    fill_vector( event, "GenLep", lep_all_varnames, leps)

    # MET
    event.GenMet_pt = reader.products['genMET'][0].pt()
    event.GenMet_phi = reader.products['genMET'][0].phi()

    # jets
    jets = map( lambda t:{var: getattr(t, var)() for var in jet_read_varnames}, filter( lambda j:j.pt()>30, reader.products['genJets']) )

    # jet/lepton disambiguation
    jets = filter( lambda j: (min([999]+[deltaR2(j, l) for l in leps if l['pt']>10]) > 0.3**2 ), jets )

    # find b's from tops:
    b_partons = [ b for b in filter( lambda p:abs(p.pdgId())==5 and p.numberOfMothers()==1 and abs(p.mother(0).pdgId())==6,  gp) ]

    for jet in jets:
        jet['matchBParton'] = ( min([999]+[deltaR2(jet, {'eta':b.eta(), 'phi':b.phi()}) for b in b_partons]) < 0.2**2 )

    jets.sort( key = lambda p:-p['pt'] )


    fill_vector( event, "GenJet", jet_write_varnames, jets)

tmp_dir     = ROOT.gDirectory
post_fix = '_%i'%args.job if args.nJobs > 1 else ''
output_filename =  os.path.join(output_directory, sample.name + post_fix + '.root')

print output_filename.replace('.root', '.log'), output_filename.replace('.root', '_rt.log')
_logger.   add_fileHandler( output_filename.replace('.root', '.log'), args.logLevel )
_logger_rt.add_fileHandler( output_filename.replace('.root', '_rt.log'), args.logLevel )

if os.path.exists( output_filename ) and not args.overwrite:
    logger.info( "File %s found. Quit.", output_filename )
    sys.exit(0)

if args.delphes:
    delphesProducer = DelphesProducer()
    delphes_file = os.path.join( output_directory, 'delphes', sample.name+'.root' )
    delphesProducer.produce( sample.files, delphes_file )

output_file = ROOT.TFile( output_filename, 'recreate')
output_file.cd()
maker = TreeMaker(
    sequence  = [ filler ],
    variables = [ TreeVariable.fromString(x) for x in variables ] + extra_variables,
    treeName = "Events"
    )

tmp_dir.cd()

counter = 0
reader.start()
maker.start()

while reader.run( ):
    #if abs(map( lambda p: p.daughter(0).pdgId(), filter( lambda p: p.pdgId()==23 and p.numberOfDaughters()==2, reader.products['gp']))[0])==13: 
    #    maker.run()
    #    break
    maker.run()

    counter += 1
    if counter == maxEvents:  break

logger.info( "Done with running over %i events.", reader.nEvents )

output_file.cd()
maker.tree.Write()
output_file.Close()

logger.info( "Written output file %s", output_filename )
