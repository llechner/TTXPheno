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
from TTXPheno.Tools.helpers                import deltaPhi, deltaR, deltaR2, cosThetaStar, closestOSDLMassToMZ, nanJet, nanLepton
from TTXPheno.Tools.HyperPoly              import HyperPoly
from TTXPheno.Tools.WeightInfo             import WeightInfo
from TTXPheno.Tools.DelphesProducer        import DelphesProducer
from TTXPheno.Tools.DelphesReader          import DelphesReader
from TTXPheno.Tools.objectSelection        import isGoodGenJet, isGoodGenLepton, isGoodRecoMuon, isGoodRecoElectron, isGoodRecoJet, isGoodRecoPhoton, genJetId

#
# Arguments
# 
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',           action='store',      default='INFO',          nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--small',              action='store_true', help='Run only on a small subset of the data?')#, default = True)
argParser.add_argument('--delphes',            action='store_true', help='Run Delphes?')
argParser.add_argument('--overwrite',          action='store_true', help='Overwrite?')#, default = True)
argParser.add_argument('--targetDir',          action='store',      default='v5')
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
    maxEvents=1000 # Number of files
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

def vecSumPt(*args):
    return sqrt( sum([o['pt']*cos(o['phi']) for o in args],0.)**2 + sum([o['pt']*sin(o['phi']) for o in args],0.)**2 )

def addIndex( collection ):
    for i  in range(len(collection)):
        collection[i]['index'] = i

# standard variables
variables  = ["run/I", "lumi/I", "evt/l"]

# MET
variables += ["genMet_pt/F", "genMet_phi/F"]

# jet vector
jet_read_vars       =  "pt/F,eta/F,phi/F,isMuon/I,isElectron/I,isPhoton/I"
jet_read_varnames   =  varnames( jet_read_vars )
jet_write_vars      = jet_read_vars+',matchBParton/I' 
jet_write_varnames  =  varnames( jet_write_vars )
variables += ["genJet[%s]"%jet_write_vars]
variables += ["genBj0_%s"%var for var in jet_write_vars.split(',')]
variables += ["genBj1_%s"%var for var in jet_write_vars.split(',')]
# lepton vector 
lep_vars       =  "pt/F,eta/F,phi/F,pdgId/I"
lep_extra_vars =  "motherPdgId/I"
lep_varnames   =  varnames( lep_vars ) 
lep_all_varnames = lep_varnames + varnames(lep_extra_vars)
variables     += ["genLep[%s]"%(','.join([lep_vars, lep_extra_vars]))]
# associated jet indices
variables += [ "genBjLeadlep_index/I", "genBjLeadhad_index/I" ]
variables += [ "genBjNonZlep_index/I", "genBjNonZhad_index/I" ]
# top vector
top_vars       =  "pt/F,eta/F,phi/F"
top_varnames   =  varnames( top_vars ) 
variables     += ["genTop[%s]"%top_vars]

# to be stored for each boson
boson_read_varnames= [ 'pt', 'phi', 'eta', 'mass']
# Z vector from gen collection
variables     += ["genZ_pt/F", "genZ_phi/F", "genZ_eta/F", "genZ_mass/F", "genZ_cosThetaStar/F", "genZ_daughterPdg/I"]
# Z vector from genleps
variables     += ["genLepZ_pt/F", "genLepZ_phi/F", "genLepZ_eta/F", "genLepZ_mass/F", "genLepZ_lldPhi/F", "genLepZ_lldR/F","genLepZ_cosThetaStar/F", "genLepZ_daughterPdg/I", "genLepNonZ_l1_index/I"]
variables     += ["genLepZ_l1_index/I", "genLepZ_l2_index/I", "genLepNonZ_l1_index/I", "genLepNonZ_l2_index/I"]
# W vector
variables     += ["genW_pt/F", "genW_phi/F", "genW_eta/F", "genW_mass/F", "genW_daughterPdg/I"]
# gamma vector
gen_photon_vars = "pt/F,phi/F,eta/F,mass/F,motherPdgId/I"
variables     += ["genPhoton[%s]"%gen_photon_vars]
gen_photon_varnames = varnames( gen_photon_vars )

if args.delphes:
    # reconstructed bosons
    variables     += ["recoZ_l1_index/I", "recoZ_l2_index/I", "recoNonZ_l1_index/I", "recoNonZ_l2_index/I",  "recoZ_pt/F", "recoZ_eta/F", "recoZ_phi/F", "recoZ_mass/F", "recoZ_lldPhi/F", "recoZ_lldR/F", "recoZ_cosThetaStar/F"]

    # reconstructed leptons
    recoLep_vars       = "pt/F,eta/F,phi/F,pdgId/I,isolationVar/F,isolationVarRhoCorr/F,sumPtCharged/F,sumPtNeutral/F,sumPtChargedPU/F,sumPt/F,ehadOverEem/F"
    variables         += ["recoLep[%s]"%recoLep_vars]
    recoLep_varnames  = varnames( recoLep_vars )
        
    # reconstructed jets
    recoJet_vars    = 'pt/F,eta/F,phi/F,bTag/F,bTagPhys/F' 
    variables      += ["recoJet[%s]"%recoJet_vars]
    recoJet_write_varnames = varnames( recoJet_vars )
    variables += ["recoBj0_%s"%var for var in recoJet_vars.split(',')]
    variables += ["recoBj1_%s"%var for var in recoJet_vars.split(',')]
    recoJet_varnames= varnames( recoJet_vars )
    
    # associated jet indices
    variables += [ "recoBjNonZlep_index/I", "recoBjNonZhad_index/I" ]
    variables += [ "recoBjLeadlep_index/I", "recoBjLeadhad_index/I" ]
    # reconstructed photons
    recoPhoton_vars = 'pt/F,eta/F,phi/F,isolationVar/F,isolationVarRhoCorr/F,sumPtCharged/F,sumPtNeutral/F,sumPtChargedPU/F,sumPt/F,ehadOverEem/F,genIndex/I'
    variables      += ["recoPhoton[%s]"%recoPhoton_vars]
    recoPhoton_varnames = varnames( recoPhoton_vars )

    variables      += ["recoMet_pt/F", "recoMet_phi/F"]
 
# Lumi weight 1fb
variables += ["lumiweight1fb/F"]
if args.addReweights:
    variables.append('rw_nominal/F')
    # Lumi weight 1fb / w_0
    variables.append("ref_lumiweight1fb/F")


def fill_vector_collection( event, collection_name, collection_varnames, objects):
    setattr( event, "n"+collection_name, len(objects) )
    for i_obj, obj in enumerate(objects):
        for var in collection_varnames:
            getattr(event, collection_name+"_"+var)[i_obj] = obj[var]
def fill_vector( event, collection_name, collection_varnames, obj):
    for var in collection_varnames:
        setattr(event, collection_name+"_"+var, obj[var] )

reader = sample.fwliteReader( products = products )

def filler( event ):

    event.run, event.lumi, event.evt = reader.evt
    event.lumiweight1fb = lumiweight1fb

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
        event.ref_lumiweight1fb = event.lumiweight1fb / coeff[0]

    # All gen particles
    gp      = reader.products['gp']

    # for searching
    search  = GenSearch( gp )

    # find heavy objects before they decay
    genTops = map( lambda t:{var: getattr(t, var)() for var in top_varnames}, filter( lambda p:abs(p.pdgId())==6 and search.isLast(p),  gp) )

    genTops.sort( key = lambda p:-p['pt'] )
    fill_vector_collection( event, "genTop", top_varnames, genTops ) 

    # generated Z's
    genZs = filter( lambda p:abs(p.pdgId())==23 and search.isLast(p) and abs(p.daughter(0).pdgId()) in [11, 13], gp)
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

    # MET
    genMet = {'pt':reader.products['genMET'][0].pt(), 'phi':reader.products['genMET'][0].phi()}
    event.genMet_pt  = genMet['pt']
    event.genMet_phi = genMet['phi'] 

    # gen photons
    genPhotons = [ (search.ascend(l), l) for l in filter( lambda p:abs(p.pdgId())==22 and p.pt()>15 and search.isLast(p), gp) ]
    genPhotons.sort( key = lambda p: -p[1].pt() )
    genPhotons_   = [] 
    for first, last in genPhotons[:100]: 
        mother_pdgId = first.mother(0).pdgId() if first.numberOfMothers()>0 else -1
        genPhotons_.append( {var: getattr(last, var)() for var in boson_read_varnames} )
        genPhotons_[-1]['motherPdgId'] = mother_pdgId

    fill_vector_collection( event, "genPhoton", gen_photon_varnames, genPhotons_ ) 
    
    # find all genLeptons 
    genLeptons = [ (search.ascend(l), l) for l in filter( lambda p:abs(p.pdgId()) in [11, 13] and search.isLast(p) and p.pt()>=0,  gp) ]
    genLeps    = []
    for first, last in genLeptons:
        mother_pdgId = first.mother(0).pdgId() if first.numberOfMothers()>0 else -1
        genLeps.append( {var: getattr(last, var)() for var in lep_varnames} )
        genLeps[-1]['motherPdgId'] = mother_pdgId

    # filter gen leptons
    genLeps =  list( filter( lambda l:isGoodGenLepton( l ), genLeps ) )
    genLeps.sort( key = lambda p:-p['pt'] )
    addIndex( genLeps )

    # jets
    fwlite_genJets = filter( genJetId, reader.products['genJets'] )
    genJets = map( lambda t:{var: getattr(t, var)() for var in jet_read_varnames}, filter( lambda j:j.pt()>30, fwlite_genJets) )
    # filter genJets
    genJets = list( filter( lambda j:isGoodGenJet( j ), genJets ) )

    # find b's from tops:
    b_partons = [ b for b in filter( lambda p:abs(p.pdgId())==5 and p.numberOfMothers()==1 and abs(p.mother(0).pdgId())==6,  gp) ]

    # store if gen-jet is DR matched to a B parton
    for genJet in genJets:
        genJet['matchBParton'] = ( min([999]+[deltaR2(genJet, {'eta':b.eta(), 'phi':b.phi()}) for b in b_partons]) < 0.2**2 )

    genJets = filter( lambda j: (min([999]+[deltaR2(j, l) for l in genLeps if l['pt']>10]) > 0.3**2 ), genJets )
    genJets.sort( key = lambda p:-p['pt'] )
    addIndex( genJets )

    # gen b jets
    trueBjets = list( filter( lambda j: j['matchBParton'], genJets ) )
    trueNonBjets = list( filter( lambda j: not j['matchBParton'], genJets ) )

    # Mimick b reconstruction ( if the trailing b fails acceptance, we supplement with the leading non-b jet ) 
    genBj0, genBj1 = ( trueBjets + trueNonBjets + [nanJet(), nanJet()] )[:2]
    if genBj0['pt']<float('inf'): fill_vector( event, "genBj0", jet_write_varnames, genBj0) 
    if genBj1['pt']<float('inf'): fill_vector( event, "genBj1", jet_write_varnames, genBj1) 

    # reco-bjet/leading lepton association
    if len(genLeps)>0 and genBj0['pt']<float('inf') and genBj1['pt']<float('inf'):
        if vecSumPt( genBj0, genLeps[0], genMet ) > vecSumPt( genBj1, genLeps[0], genMet ):
            event.genBjLeadlep_index, event.genBjLeadhad_index = genBj0['index'], genBj1['index']
        else:
            event.genBjLeadlep_index, event.genBjLeadhad_index = genBj1['index'], genBj0['index']

    # find Z in genLep
    (event.genLepZ_mass, genLepZ_l1_index, genLepZ_l2_index) = closestOSDLMassToMZ(genLeps)
    genLepNonZ_indices = [ i for i in range(len(genLeps)) if i not in [genLepZ_l1_index, genLepZ_l2_index] ]
    event.genLepZ_l1_index    = genLeps[genLepZ_l1_index]['index'] if genLepZ_l1_index>=0 else -1
    event.genLepZ_l2_index    = genLeps[genLepZ_l2_index]['index'] if genLepZ_l2_index>=0 else -1
    event.genLepNonZ_l1_index = genLeps[genLepNonZ_indices[0]]['index'] if len(genLepNonZ_indices)>0 else -1
    event.genLepNonZ_l2_index = genLeps[genLepNonZ_indices[1]]['index'] if len(genLepNonZ_indices)>1 else -1
    # store genLepZ stuff
    if event.genLepZ_mass>0:
        genLepZ_l1 = ROOT.TLorentzVector()
        genLepZ_l1.SetPtEtaPhiM(genLeps[event.genLepZ_l1_index]['pt'], genLeps[event.genLepZ_l1_index]['eta'], genLeps[event.genLepZ_l1_index]['phi'], 0 )
        genLepZ_l2 = ROOT.TLorentzVector()
        genLepZ_l2.SetPtEtaPhiM(genLeps[event.genLepZ_l2_index]['pt'], genLeps[event.genLepZ_l2_index]['eta'], genLeps[event.genLepZ_l2_index]['phi'], 0 )
        genLepZ = genLepZ_l1 + genLepZ_l2
        event.genLepZ_pt   = genLepZ.Pt()
        event.genLepZ_eta  = genLepZ.Eta()
        event.genLepZ_phi  = genLepZ.Phi()
        event.genLepZ_lldPhi = deltaPhi(genLeps[event.genLepZ_l1_index]['phi'], genLeps[event.genLepZ_l2_index]['phi'])
        event.genLepZ_lldR   = deltaR(genLeps[event.genLepZ_l1_index], genLeps[event.genLepZ_l2_index])
        genLepMinus_index = event.genLepZ_l1_index if genLeps[event.genLepZ_l1_index]['pdgId'] > 0 else event.genLepZ_l2_index
        event.genLepZ_cosThetaStar = cosThetaStar(event.genLepZ_mass, event.genLepZ_pt, event.genLepZ_eta, event.genLepZ_phi, genLeps[genLepMinus_index]['pt'], genLeps[genLepMinus_index]['eta'], genLeps[genLepMinus_index]['phi'] )

    # reco-bjet/nonZ lepton association
    if event.genLepNonZ_l1_index>=0 and genBj0['pt']<float('inf') and genBj1['pt']<float('inf'):
        if vecSumPt( genBj0, genLeps[event.genLepNonZ_l1_index], genMet ) > vecSumPt( genBj1, genLeps[event.genLepNonZ_l1_index], genMet ):
            event.genBjNonZlep_index, event.genBjNonZhad_index = genBj0['index'], genBj1['index']
        else:
            event.genBjNonZlep_index, event.genBjNonZhad_index = genBj1['index'], genBj0['index']

    #for jet in genJets:
    #    print jet['isMuon'], jet['isElectron'], jet['isPhoton'], min([999]+[deltaR2(jet, l) for l in genLeps if l['pt']>10]), jet

    # jet/lepton disambiguation -> remove jets, because gen-jets cluster all leptons
    #if args.logLevel == 'DEBUG':
    #    for jet in filter( lambda j: not (min([999]+[deltaR2(j, l) for l in genLeps if l['pt']>10]) > 0.3**2 ), genJets ):
    #        logger.debug( "Filtered gen %f jet %r lep %r", sqrt((min([999]+[deltaR2(jet, l) for l in genLeps if l['pt']>10]))), jet, [ (l['eta'], jet['pt']/l['pt']) for l in genLeps] )
    #        assert False, ""


    fill_vector_collection( event, "genLep", lep_all_varnames, genLeps)
    fill_vector_collection( event, "genJet", jet_write_varnames, genJets)

    # Reco quantities
    if args.delphes:
        delphesReader.getEntry(reader.position-1 )

        # read jets
        recoJets =  filter( isGoodRecoJet, delphesReader.jets()) 
        recoJets.sort( key = lambda p:-p['pt'] )
        addIndex( recoJets )

        # make reco b jets
        recoBJets    = filter( lambda j:j['bTag']==1, recoJets )
        recoNonBJets = filter( lambda j:not (j['bTag']==1), recoJets )
        recoBj0, recoBj1 = ( recoBJets + recoNonBJets + [nanJet(), nanJet()] )[:2] 
        fill_vector( event, "recoBj0", recoJet_write_varnames, recoBj0) 
        fill_vector( event, "recoBj1", recoJet_write_varnames, recoBj1) 

        # read leptons
        recoLeps =  filter( isGoodRecoMuon, delphesReader.muons()) + filter( isGoodRecoElectron, delphesReader.electrons() )
        recoLeps.sort( key = lambda p:-p['pt'] )

        for i_lep, lep in enumerate(recoLeps):
            min_DR=999
            min_jet=None
            for jet in recoJets:
                if deltaR2(jet, lep)<min_DR:
                    min_DR_jet = jet
                    min_DR     = deltaR2(jet, lep)
            if min_DR<0.4**2: 
                #print "Filtering lepton", lep, "because DR", sqrt(min_DR), "of jet", jet
                break

        # cross-cleaning of reco-objects
        recoLeps = filter( lambda l: (min([999]+[deltaR2(l, j) for j in recoJets if j['pt']>30]) > 0.3**2 ), recoLeps )
        # give index to leptons
        addIndex( recoLeps )

        # MET
        recoMet = delphesReader.met()[0]

        # reco-bjet/leading lepton association
        if len(recoLeps)>0 and recoBj0['pt']<float('inf') and recoBj1['pt']<float('inf'):
            if vecSumPt( recoBj0, recoLeps[0], recoMet ) > vecSumPt( recoBj1, recoLeps[0], recoMet ):
                event.recoBjLeadlep_index, event.recoBjLeadhad_index = recoBj0['index'], recoBj1['index']
            else:
                event.recoBjLeadlep_index, event.recoBjLeadhad_index = recoBj1['index'], recoBj0['index']

        # Photons
        recoPhotons = filter( isGoodRecoPhoton, delphesReader.photons() )
        for recoPhoton in recoPhotons:
            recoPhoton['genIndex'] = -1
            minDR = 999
            for index, genPhoton in enumerate(genPhotons_):
                dr = deltaR( recoPhoton, genPhoton )
                if dr<0.4 and dr<minDR:
                    minDR = dr
                    recoPhoton['genIndex'] = index 
        # Store
        fill_vector_collection( event, "recoLep",    recoLep_varnames, recoLeps )
        fill_vector_collection( event, "recoJet",    recoJet_varnames, recoJets )
        fill_vector_collection( event, "recoPhoton", recoPhoton_varnames, recoPhotons )

        event.recoMet_pt  = recoMet['pt']
        event.recoMet_phi = recoMet['phi']

        # search for reco Z in reco leptons
        (event.recoZ_mass, recoZ_l1_index, recoZ_l2_index) = closestOSDLMassToMZ(recoLeps)
        recoNonZ_indices = [ i for i in range(len(recoLeps)) if i not in [recoZ_l1_index, recoZ_l2_index] ]
        event.recoZ_l1_index    = recoLeps[recoZ_l1_index]['index'] if recoZ_l1_index>=0 else -1
        event.recoZ_l2_index    = recoLeps[recoZ_l2_index]['index'] if recoZ_l2_index>=0 else -1
        event.recoNonZ_l1_index = recoLeps[recoNonZ_indices[0]]['index'] if len(recoNonZ_indices)>0 else -1
        event.recoNonZ_l2_index = recoLeps[recoNonZ_indices[1]]['index'] if len(recoNonZ_indices)>1 else -1

        # Store Z information 
        if event.recoZ_mass>=0:
            if recoLeps[event.recoZ_l1_index]['pdgId']*recoLeps[event.recoZ_l2_index]['pdgId']>0 or abs(recoLeps[event.recoZ_l1_index]['pdgId'])!=abs(recoLeps[event.recoZ_l2_index]['pdgId']): 
                raise RuntimeError( "not a Z! Should not happen" )
            Z_l1 = ROOT.TLorentzVector()
            Z_l1.SetPtEtaPhiM(recoLeps[event.recoZ_l1_index]['pt'], recoLeps[event.recoZ_l1_index]['eta'], recoLeps[event.recoZ_l1_index]['phi'], 0 )
            Z_l2 = ROOT.TLorentzVector()
            Z_l2.SetPtEtaPhiM(recoLeps[event.recoZ_l2_index]['pt'], recoLeps[event.recoZ_l2_index]['eta'], recoLeps[event.recoZ_l2_index]['phi'], 0 )
            Z = Z_l1 + Z_l2
            event.recoZ_pt   = Z.Pt()
            event.recoZ_eta  = Z.Eta()
            event.recoZ_phi  = Z.Phi()
            event.recoZ_lldPhi = deltaPhi(recoLeps[event.recoZ_l1_index]['phi'], recoLeps[event.recoZ_l2_index]['phi'])
            event.recoZ_lldR   = deltaR(recoLeps[event.recoZ_l1_index], recoLeps[event.recoZ_l2_index])
            lm_index = event.recoZ_l1_index if recoLeps[event.recoZ_l1_index]['pdgId'] > 0 else event.recoZ_l2_index
            event.recoZ_cosThetaStar = cosThetaStar(event.recoZ_mass, event.recoZ_pt, event.recoZ_eta, event.recoZ_phi, recoLeps[lm_index]['pt'], recoLeps[lm_index]['eta'], recoLeps[lm_index]['phi'] )

            # reco-bjet/lepton association
            if event.recoNonZ_l1_index>=0 and recoBj0['pt']<float('inf') and recoBj1['pt']<float('inf'):
                if vecSumPt( recoBj0, recoLeps[event.recoNonZ_l1_index], recoMet ) > vecSumPt( recoBj1, recoLeps[event.recoNonZ_l1_index], recoMet ):
                    event.recoBjNonZlep_index, event.recoBjNonZhad_index = recoBj0['index'], recoBj1['index']
                else:
                    event.recoBjNonZlep_index, event.recoBjNonZhad_index = recoBj1['index'], recoBj0['index']

         
tmp_dir     = ROOT.gDirectory
#post_fix = '_%i'%args.job if args.nJobs > 1 else ''
output_filename =  os.path.join(output_directory, sample.name + '.root')

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
