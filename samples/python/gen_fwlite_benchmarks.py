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
