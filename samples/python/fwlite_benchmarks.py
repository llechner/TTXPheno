''' Benchmark samples for TTXPheno (EDM)'''

# standard imports
import os

# RootTools
from RootTools.core.standard import *

# TTXPheno
from TTXPheno.Tools.user import results_directory 

# sqlite3 sample cache file
dbFile = os.path.join( results_directory, 'sample_cache', 'fwlite_benchmarks.db')

# Logging
if __name__ == "__main__":
    import TTXPheno.Tools.logger as logger
    logger = logger.get_logger('DEBUG')
    import RootTools.core.logger as logger_rt
    logger_rt = logger_rt.get_logger('DEBUG')

# TTZ 931 gen with weights - high stats!
fwlite_ttZ_ll_LO_highStat_scan     = FWLiteSample.fromDAS("fwlite_ttZ_ll_LO_highStat_scan", "/ttZ0j_rwgt_patch_625_highStat_slc6_amd64_gcc630_CMSSW_9_3_0_tarball/schoef-dim6top_07April18-37d73f7f997f18e72dbfd34806877f87/USER", "phys03", dbFile = dbFile)
fwlite_ttZ_ll_LO_highStat_scan.reweight_pkl = "/afs/hephy.at/data/rschoefbeck02/TopEFT/results/gridpacks/ttZ0j_rwgt_patch_625_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"

# TTZ 931 gen with weights - current plane high stats!
fwlite_ttZ_ll_LO_currentplane_highStat_scan     = FWLiteSample.fromDAS("fwlite_ttZ_ll_LO_currentplane_highStat_scan", "/ttZ0j_rwgt_patch_currentplane_highStat_slc6_amd64_gcc630_CMSSW_9_3_0_tarball/schoef-dim6top_07April18-b9d79c4e9bec84f06f452ff39977163a/USER", "phys03", dbFile = dbFile)
fwlite_ttZ_ll_LO_currentplane_highStat_scan.reweight_pkl = "/afs/hephy.at/data/rschoefbeck02/TopEFT/results/gridpacks/ttZ0j_rwgt_patch_currentplane_highStat_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"

# test for Lukas
test = FWLiteSample.fromFiles("test", files = ["/afs/hephy.at/work/r/rschoefbeck/CMS/gen/CMSSW_9_3_1/src/TopEFT/Generation/production/GEN_LO_0j.root"])
test.reweight_pkl = '/afs/hephy.at/data/llechner01/TTXPheno/gridpacks/07052018/ttZ/order3/ttZ0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl'
