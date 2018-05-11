''' Benchmark samples for TopEFT (EDM)'''

# standard imports
import os

# RootTools
from RootTools.core.standard import *

#Top EFT
from TTXPheno.Tools.user import results_directory 

# Logging
import logging
logger = logging.getLogger(__name__)

import glob

gen_dir = "/afs/hephy.at/data/rschoefbeck02/TopEFT/skims/gen/v2/"

# Robert first ttZ_ll dim6top scan
dim6top_ttZ_ll_LO_highStat_scan                           = Sample.fromDirectory("fwlite_ttZ_ll_LO_highStat_scan", texName = "ttZ (scan)", directory = [os.path.join( gen_dir, "fwlite_ttZ_ll_LO_highStat_scan")])
dim6top_ttZ_ll_LO_highStat_scan.reweight_pkl              = "/afs/hephy.at/data/rschoefbeck02/TopEFT/results/gridpacks/ttZ0j_rwgt_patch_625_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"

dim6top_ttZ_ll_LO_currentplane_highStat_scan              = Sample.fromDirectory("dim6top_ttZ_ll_LO_currentplane_highStat_scan", texName = "ttZ (current scan)", directory = [os.path.join( gen_dir, "fwlite_ttZ_ll_LO_currentplane_highStat_scan")])
dim6top_ttZ_ll_LO_currentplane_highStat_scan.reweight_pkl = "/afs/hephy.at/data/rschoefbeck02/TopEFT/results/gridpacks/ttZ0j_rwgt_patch_currentplane_highStat_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"

test = Sample.fromFiles("test", files = ["/afs/hephy.at/data/rschoefbeck02/TTXPheno/skims/gen/v2/test/test.root"], texName = "test")
test.reweight_pkl = '/afs/cern.ch/user/l/llechner/public/gridpacks_data/order_3/ttZ/gridpacks/reweight_card.pkl'

gen_dir = "/afs/hephy.at/data/rschoefbeck01/TTXPheno/skims/gen/v3/"
gridpack_dir = "/afs/hephy.at/data/llechner01/TTXPheno/gridpacks/"

fwlite_ttZ_ll_LO_order3_8weights  = Sample.fromDirectory("fwlite_ttZ_ll_LO_order3_8weights",  texName = "ttZ",      directory = [os.path.join( gen_dir, "fwlite_ttZ_ll_LO_order3_8weights" )]) 
fwlite_ttZ_ll_LO_order3_8weights.reweight_pkl = gridpack_dir + "07052018/ttZ/order3/ttZ0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"

fwlite_ttW_LO_order3_8weights     = Sample.fromDirectory("fwlite_ttW_LO_order3_8weights",     texName = "ttW",      directory = [os.path.join( gen_dir, "fwlite_ttW_LO_order3_8weights" )])
fwlite_ttW_ll_LO_order3_8weights.reweight_pkl = gridpack_dir + "07052018/ttW/order3/ttW0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"

fwlite_ttgamma_LO_order3_8weights = Sample.fromDirectory("fwlite_ttgamma_LO_order3_8weights", texName = "ttgamma",  directory = [os.path.join( gen_dir, "fwlite_ttgamma_LO_order3_8weights" )])
fwlite_ttgamma_ll_LO_order3_8weights.reweight_pkl = gridpack_dir + "07052018/ttgamma/order3/ttgamma0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"

fwlite_ttZ_ll_LO_order3_8weights_delphes  = Sample.fromDirectory("fwlite_ttZ_ll_LO_order3_8weights_delphes",  texName = "ttZ",      directory = [os.path.join( gen_dir, "fwlite_ttZ_ll_LO_order3_8weights", "delphes" )], treeName = "Delphes") 
fwlite_ttW_LO_order3_8weights_delphes     = Sample.fromDirectory("fwlite_ttW_LO_order3_8weights_delphes",     texName = "ttW",      directory = [os.path.join( gen_dir, "fwlite_ttW_LO_order3_8weights", "delphes" )], treeName = "Delphes")
fwlite_ttgamma_LO_order3_8weights_delphes = Sample.fromDirectory("fwlite_ttgamma_LO_order3_8weights_delphes", texName = "ttgamma",  directory = [os.path.join( gen_dir, "fwlite_ttgamma_LO_order3_8weights", "delphes" )], treeName = "Delphes")

fwlite_ttZ_ll_LO_order3_8weights.addFriend(  fwlite_ttZ_ll_LO_order3_8weights_delphes, treeName = "Delphes", sortFiles = True)
fwlite_ttW_LO_order3_8weights.addFriend(  fwlite_ttW_LO_order3_8weights_delphes, treeName = "Delphes", sortFiles = True)
fwlite_ttgamma_LO_order3_8weights.addFriend(  fwlite_ttgamma_LO_order3_8weights_delphes, treeName = "Delphes", sortFiles = True)
