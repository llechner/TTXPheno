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
from TTXPheno.Tools.DelphesReader          import DelphesReader

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

maxEvents = -1
if args.small: 
    args.targetDir += "_small"
    maxEvents=10 # Number of files
    sample.files=sample.files[:1]

xsec = sample.xsec
nEvents = sample.nEvents
lumiweight1fb = xsec * 1000. / nEvents

# output directory
output_directory = os.path.join(skim_output_directory, 'gen', args.targetDir, sample.name) 
if not os.path.exists( output_directory ): 
    os.makedirs( output_directory )
    logger.info( "Created output directory %s", output_directory )

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
variables += ["genMet_pt/F", "genMet_phi/F"]

# jet vector
jet_read_vars       =  "pt/F,eta/F,phi/F"
jet_read_varnames   =  varnames( jet_read_vars )
jet_write_vars      = jet_read_vars+',matchBParton/I' 
jet_write_varnames  =  varnames( jet_write_vars )
variables     += ["genJet[%s]"%jet_write_vars]
# lepton vector 
lep_vars       =  "pt/F,eta/F,phi/F,pdgId/I"
lep_extra_vars =  "motherPdgId/I"
lep_varnames   =  varnames( lep_vars ) 
lep_all_varnames = lep_varnames + varnames(lep_extra_vars)
variables     += ["genLep[%s]"%(','.join([lep_vars, lep_extra_vars]))]
# top vector
top_vars       =  "pt/F,eta/F,phi/F"
top_varnames   =  varnames( top_vars ) 
variables     += ["genTop[%s]"%top_vars]

# to be stored for each boson
boson_read_varnames= [ 'pt', 'phi', 'eta', 'mass']
# Z vector
variables     += ["genZ_pt/F", "genZ_phi/F", "genZ_eta/F", "genZ_mass/F", "genZ_cosThetaStar/F", "genZ_daughterPdg/I"]
# W vector
variables     += ["genW_pt/F", "genW_phi/F", "genW_eta/F", "genW_mass/F", "genW_daughterPdg/I"]
# gamma vector
variables     += ["genGamma_pt/F", "genGamma_phi/F", "genGamma_eta/F", "genGamma_mass/F"]

# Lumi weight 1fb
variables += ["lumiweight1fb/F"]
if args.addReweights:
    variables.append('rw_nominal/F')
    # Lumi weight 1fb / w_0
    variables.append("ref_lumiweight1fb/F")


def fill_vector( event, collection_name, collection_varnames, objects):
    setattr( event, "n"+collection_name, len(objects) )
    for i_obj, obj in enumerate(objects):
        for var in collection_varnames:
            getattr(event, collection_name+"_"+var)[i_obj] = obj[var]

reader = sample.fwliteReader( products = products )

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

        # get list of values of ref point in specific order
        ref_point_coordinates = [weightInfo.ref_point_coordinates[var] for var in weightInfo.variables]

        # Initialize with Reference Point
        if not hyperPoly.initialized: hyperPoly.initialize( param_points, ref_point_coordinates )
        coeff = hyperPoly.get_parametrization( weights )

        # = HyperPoly(weight_data, args.interpolationOrder)
        event.np = hyperPoly.ndof
        event.chi2_ndof = hyperPoly.chi2_ndof(coeff, weights)
        #logger.debug( "chi2_ndof %f coeff %r", event.chi2_ndof, coeff )
        if event.chi2_ndof>10**-6: logger.warning( "chi2_ndof is large: %f", event.chi2_ndof )
        for n in xrange(hyperPoly.ndof):
            event.p_C[n] = coeff[n]

        # lumi weight / w0
        event.lumiweight1fb = lumiweight1fb
        if args.addReweights:
            event.ref_lumiweight1fb = event.lumiweight1fb / coeff[0]

    # All gen particles
    gp      = reader.products['gp']

    # for searching
    search  = GenSearch( gp )

    # find heavy objects before they decay
    genTops = map( lambda t:{var: getattr(t, var)() for var in top_varnames}, filter( lambda p:abs(p.pdgId())==6 and search.isLast(p),  gp) )

    genTops.sort( key = lambda p:-p['pt'] )
    fill_vector( event, "genTop", top_varnames, genTops ) 

    # generated Z's
    genZs = filter( lambda p:abs(p.pdgId())==23 and search.isLast(p), gp)
    genZs.sort( key = lambda p: -p.pt() )
    if len(genZs)>0: 
        genZ = genZs[0]
        for var in boson_read_varnames:
           setattr( event, "genZ_"+var,  getattr(genZ, var)() )
    else:
        genZ = None
    
    if genZ is not None:

        d1, d2 = genZ.daughter(0), genZ.daughter(1)
        if d1.pdgId()>0: 
            lm, lp = d1, d2
        else:
            lm, lp = d2, d1
        event.genZ_daughterPdg = lm.pdgId()
        event.genZ_cosThetaStar = cosThetaStar(genZ.mass(), genZ.pt(), genZ.eta(), genZ.phi(), lm.pt(), lm.eta(), lm.phi())

    # generated W's
    genWs = filter( lambda p:abs(p.pdgId())==24 and search.isLast(p), gp)
    genWs.sort( key = lambda p: -p.pt() )
    # W can't have a top-mother - We're looking for the extra boson (there is always an extra boson)
    genWs = filter( lambda p: abs(search.ascend(p).mother(0).pdgId())!=6, genWs )
    if len(genWs)>0: 
        genW = genWs[0]
        for var in boson_read_varnames:
           setattr( event, "genW_"+var,  getattr(genW, var)() )
    else:
        genW = None
    
    if genW is not None:

        d1, d2 = genW.daughter(0), genW.daughter(1)
        if abs(d1.pdgId()) in [11, 13, 15]: 
            lep, neu = d1, d2
        else:
            lep, neu = d2, d1

        event.genW_daughterPdg = lep.pdgId()

    genGammas = filter( lambda p:abs(p.pdgId())==22 and search.isLast(p), gp)
    genGammas.sort( key = lambda p: -p.pt() )
    if len(genGammas)>0: 
        genGamma = genGammas[0]
        for var in boson_read_varnames:
           setattr( event, "genGamma_"+var,  getattr(genGamma, var)() )
    else:
        genGamma = None
    
    # find all genLeptons 
    genLeptons = [ (search.ascend(l), l) for l in filter( lambda p:abs(p.pdgId()) in [11, 13] and search.isLast(p) and p.pt()>=0,  gp) ]
    genLeps    = []
    for first, last in genLeptons:
        mother_pdgId = first.mother(0).pdgId() if first.numberOfMothers()>0 else -1
        genLeps.append( {var: getattr(last, var)() for var in lep_varnames} )
        genLeps[-1]['motherPdgId'] = mother_pdgId

    genLeps.sort( key = lambda p:-p['pt'] )
    fill_vector( event, "genLep", lep_all_varnames, genLeps)

    # MET
    event.genMet_pt = reader.products['genMET'][0].pt()
    event.genMet_phi = reader.products['genMET'][0].phi()

    # jets
    genJets = map( lambda t:{var: getattr(t, var)() for var in jet_read_varnames}, filter( lambda j:j.pt()>30, reader.products['genJets']) )

    # jet/lepton disambiguation
    genJets = filter( lambda j: (min([999]+[deltaR2(j, l) for l in genLeps if l['pt']>10]) > 0.3**2 ), genJets )

    # find b's from tops:
    b_partons = [ b for b in filter( lambda p:abs(p.pdgId())==5 and p.numberOfMothers()==1 and abs(p.mother(0).pdgId())==6,  gp) ]

    for genJet in genJets:
        genJet['matchBParton'] = ( min([999]+[deltaR2(genJet, {'eta':b.eta(), 'phi':b.phi()}) for b in b_partons]) < 0.2**2 )

    genJets.sort( key = lambda p:-p['pt'] )
    fill_vector( event, "genJet", jet_write_varnames, genJets)

    # Reco quantities
    #if args.delphes:
    #    delphesReader.getEntry(reader.position-1 )
    #    print delphesReader.muons()
    #    print delphesReader.electrons()
    #    print delphesReader.photons()
    #    print delphesReader.met()

tmp_dir     = ROOT.gDirectory
#post_fix = '_%i'%args.job if args.nJobs > 1 else ''
output_filename =  os.path.join(output_directory, sample.name + '.root')

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
    delphesReader = DelphesReader( delphes_file )

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
