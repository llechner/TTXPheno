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

# Logging
if __name__ == "__main__":
    import TTXPheno.Tools.logger as logger
    logger = logger.get_logger('DEBUG')
    import RootTools.core.logger as logger_rt
    logger_rt = logger_rt.get_logger('DEBUG')

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

# no reference point 8/3
fwlite_ttZ_ll_LO_order3_8weights               = Sample.fromDirectory("fwlite_ttZ_ll_LO_order3_8weights",  texName = "ttZ",      directory = [os.path.join( gen_dir, "fwlite_ttZ_ll_LO_order3_8weights" )]) 
fwlite_ttZ_ll_LO_order3_8weights.reweight_pkl  = gridpack_dir + "07052018/ttZ/order3/ttZ0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttZ_ll_LO_order3_8weights.nEvents       = 995000
fwlite_ttZ_ll_LO_order3_8weights.xsec          = 0.0915 #pb ttZ, Z->ll

fwlite_ttW_LO_order3_8weights                  = Sample.fromDirectory("fwlite_ttW_LO_order3_8weights",     texName = "ttW",      directory = [os.path.join( gen_dir, "fwlite_ttW_LO_order3_8weights" )])
fwlite_ttW_LO_order3_8weights.reweight_pkl     = gridpack_dir + "07052018/ttW/order3/ttW0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttW_LO_order3_8weights.nEvents          = 1000000
fwlite_ttW_LO_order3_8weights.xsec             = 0.2043 #pb ttW, W->lnu

fwlite_ttgamma_LO_order3_8weights              = Sample.fromDirectory("fwlite_ttgamma_LO_order3_8weights", texName = "ttgamma",  directory = [os.path.join( gen_dir, "fwlite_ttgamma_LO_order3_8weights" )])
fwlite_ttgamma_LO_order3_8weights.reweight_pkl = gridpack_dir + "07052018/ttgamma/order3/ttgamma0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttgamma_LO_order3_8weights.nEvents      = 1000000
fwlite_ttgamma_LO_order3_8weights.xsec         = 3.697 #pb ttgamma

fwlite_ttZ_ll_LO_order3_8weights_delphes       = Sample.fromDirectory("fwlite_ttZ_ll_LO_order3_8weights_delphes",  texName = "ttZ",      directory = [os.path.join( gen_dir, "fwlite_ttZ_ll_LO_order3_8weights", "delphes" )], treeName = "Delphes") 
fwlite_ttW_LO_order3_8weights_delphes          = Sample.fromDirectory("fwlite_ttW_LO_order3_8weights_delphes",     texName = "ttW",      directory = [os.path.join( gen_dir, "fwlite_ttW_LO_order3_8weights", "delphes" )], treeName = "Delphes")
fwlite_ttgamma_LO_order3_8weights_delphes      = Sample.fromDirectory("fwlite_ttgamma_LO_order3_8weights_delphes", texName = "ttgamma",  directory = [os.path.join( gen_dir, "fwlite_ttgamma_LO_order3_8weights", "delphes" )], treeName = "Delphes")

#fwlite_ttZ_ll_LO_order3_8weights.addFriend(  fwlite_ttZ_ll_LO_order3_8weights_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_ttW_LO_order3_8weights.addFriend(  fwlite_ttW_LO_order3_8weights_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_ttgamma_LO_order3_8weights.addFriend(  fwlite_ttgamma_LO_order3_8weights_delphes, treeName = "Delphes", sortFiles = True)

gen_dir = "/afs/hephy.at/data/llechner01/TTXPheno/skims/gen/v3/"

# no reference point 15/2
fwlite_ttZ_ll_LO_order2_15weights               = Sample.fromDirectory("fwlite_ttZ_ll_LO_order2_15weights",  texName = "ttZ",      directory = [os.path.join( gen_dir, "fwlite_ttZ_ll_LO_order2_15weights" )]) 
fwlite_ttZ_ll_LO_order2_15weights.reweight_pkl  = gridpack_dir + "18052018/ttZ/order2/ttZ0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttZ_ll_LO_order2_15weights.nEvents       = 970000
fwlite_ttZ_ll_LO_order2_15weights.xsec          = 0.0915 #pb ttZ, Z->ll

fwlite_ttW_LO_order2_15weights                  = Sample.fromDirectory("fwlite_ttW_LO_order2_15weights",     texName = "ttW",      directory = [os.path.join( gen_dir, "fwlite_ttW_LO_order2_15weights" )])
fwlite_ttW_LO_order2_15weights.reweight_pkl     = gridpack_dir + "18052018/ttW/order2/ttW0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttW_LO_order2_15weights.nEvents          = 1000000
fwlite_ttW_LO_order2_15weights.xsec             = 0.2043 #pb ttW, W->lnu

fwlite_ttgamma_LO_order2_15weights              = Sample.fromDirectory("fwlite_ttgamma_LO_order2_15weights", texName = "ttgamma",  directory = [os.path.join( gen_dir, "fwlite_ttgamma_LO_order2_15weights" )])
fwlite_ttgamma_LO_order2_15weights.reweight_pkl = gridpack_dir + "18052018/ttgamma/order2/ttgamma0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttgamma_LO_order2_15weights.nEvents      = 990000
fwlite_ttgamma_LO_order2_15weights.xsec         = 3.697 #pb ttgamma

fwlite_ttZ_ll_LO_order2_15weights_delphes       = Sample.fromDirectory("fwlite_ttZ_ll_LO_order2_15weights_delphes",  texName = "ttZ",      directory = [os.path.join( gen_dir, "fwlite_ttZ_ll_LO_order2_15weights", "delphes" )], treeName = "Delphes") 
fwlite_ttW_LO_order2_15weights_delphes          = Sample.fromDirectory("fwlite_ttW_LO_order2_15weights_delphes",     texName = "ttW",      directory = [os.path.join( gen_dir, "fwlite_ttW_LO_order2_15weights", "delphes" )], treeName = "Delphes")
fwlite_ttgamma_LO_order2_15weights_delphes      = Sample.fromDirectory("fwlite_ttgamma_LO_order2_15weights_delphes", texName = "ttgamma",  directory = [os.path.join( gen_dir, "fwlite_ttgamma_LO_order2_15weights", "delphes" )], treeName = "Delphes")

#fwlite_ttZ_ll_LO_order2_15weights.addFriend(  fwlite_ttZ_ll_LO_order2_15weights_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_ttW_LO_order2_15weights.addFriend(  fwlite_ttW_LO_order2_15weights_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_ttgamma_LO_order2_15weights.addFriend(  fwlite_ttgamma_LO_order2_15weights_delphes, treeName = "Delphes", sortFiles = True)

# reference point 15/2
fwlite_ttZ_ll_LO_order2_15weights_ref               = Sample.fromDirectory("fwlite_ttZ_ll_LO_order2_15weights_ref",  texName = "ttZ",      directory = [os.path.join( gen_dir, "fwlite_ttZ_ll_LO_order2_15weights_ref" )]) 
fwlite_ttZ_ll_LO_order2_15weights_ref.reweight_pkl  = gridpack_dir + "18052018_ref/ttZ/order2/ttZ0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttZ_ll_LO_order2_15weights_ref.nEvents       = 990000
fwlite_ttZ_ll_LO_order2_15weights_ref.xsec          = 0.0915 * 0.6823 / 0.07082 #pb ttZ, Z->ll, xsec_SM_NNLO * xsec_BSM_LO / xsec_SM_LO

fwlite_ttW_LO_order2_15weights_ref                  = Sample.fromDirectory("fwlite_ttW_LO_order2_15weights_ref",     texName = "ttW",      directory = [os.path.join( gen_dir, "fwlite_ttW_LO_order2_15weights_ref" )])
fwlite_ttW_LO_order2_15weights_ref.reweight_pkl     = gridpack_dir + "18052018_ref/ttW/order2/ttW0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttW_LO_order2_15weights_ref.nEvents          = 945000
fwlite_ttW_LO_order2_15weights_ref.xsec             = 0.2043 * 0.4097 / 0.1323 #pb ttW, W->lnu, xsec_SM_NNLO * xsec_BSM_LO / xsec_SM_LO

fwlite_ttgamma_LO_order2_15weights_ref              = Sample.fromDirectory("fwlite_ttgamma_LO_order2_15weights_ref", texName = "ttgamma",  directory = [os.path.join( gen_dir, "fwlite_ttgamma_LO_order2_15weights_ref" )])
fwlite_ttgamma_LO_order2_15weights_ref.reweight_pkl = gridpack_dir + "18052018_ref/ttgamma/order2/ttgamma0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttgamma_LO_order2_15weights_ref.nEvents      = 970000
fwlite_ttgamma_LO_order2_15weights_ref.xsec         = 3.697 * 7.838 / 2.439 #pb ttgamma, xsec_SM_NNLO * xsec_BSM_LO / xsec_SM_LO

fwlite_ttZ_ll_LO_order2_15weights_ref_delphes       = Sample.fromDirectory("fwlite_ttZ_ll_LO_order2_15weights_ref_delphes",  texName = "ttZ",      directory = [os.path.join( gen_dir, "fwlite_ttZ_ll_LO_order2_15weights_ref", "delphes" )], treeName = "Delphes") 
fwlite_ttW_LO_order2_15weights_ref_delphes          = Sample.fromDirectory("fwlite_ttW_LO_order2_15weights_ref_delphes",     texName = "ttW",      directory = [os.path.join( gen_dir, "fwlite_ttW_LO_order2_15weights_ref", "delphes" )], treeName = "Delphes")
fwlite_ttgamma_LO_order2_15weights_ref_delphes      = Sample.fromDirectory("fwlite_ttgamma_LO_order2_15weights_ref_delphes", texName = "ttgamma",  directory = [os.path.join( gen_dir, "fwlite_ttgamma_LO_order2_15weights_ref", "delphes" )], treeName = "Delphes")

#fwlite_ttZ_ll_LO_order2_15weights_ref.addFriend(  fwlite_ttZ_ll_LO_order2_15weights_ref_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_ttW_LO_order2_15weights_ref.addFriend(  fwlite_ttW_LO_order2_15weights_ref_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_ttgamma_LO_order2_15weights_ref.addFriend(  fwlite_ttgamma_LO_order2_15weights_ref_delphes, treeName = "Delphes", sortFiles = True)

"""
# backgrounds 15/2 (tt) or 0/2 (WZ)
# leptonic decays W > lnu, Z > ll, t > Wb
fwlite_tt_lep_LO_order2_15weights_ref               = Sample.fromDirectory("fwlite_tt_lep_LO_order2_15weights_ref",  texName = "tt",      directory = [os.path.join( gen_dir, "fwlite_tt_lep_LO_order2_15weights_ref" )]) 
fwlite_tt_lep_LO_order2_15weights_ref.reweight_pkl  = gridpack_dir + "05062018/tt/order2/tt_lep_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_tt_lep_LO_order2_15weights_ref.nEvents       = 1000000 #? not checked!
fwlite_tt_lep_LO_order2_15weights_ref.xsec          = 87.3148 #pb tt, W->lnu, xsec_SM_NNLO
fwlite_tt_lep_LO_order2_15weights_ref.xsecSM_LO     = 1 #pb tt, W->lnu not checked!
fwlite_tt_lep_LO_order2_15weights_ref.xsecBSM_LO    = 1 #pb tt, W->lnu not checked!

fwlite_WZ_lep_LO_order2_15weights                   = Sample.fromDirectory("fwlite_WZ_lep_LO_order2_15weights",     texName = "WZ",      directory = [os.path.join( gen_dir, "fwlite_WZ_lep_LO_order2_15weights" )])
fwlite_WZ_lep_LO_order2_15weights.reweight_pkl      = gridpack_dir + "05062018/WZ/order2/WZ_lep_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_WZ_lep_LO_order2_15weights.nEvents           = 1000000 #? not checked!
fwlite_WZ_lep_LO_order2_15weights.xsec              = 4.666 #pb WZ, W->lnu, Z->ll, xsec_SM_NNLO

fwlite_tt_lep_LO_order2_15weights_ref_delphes       = Sample.fromDirectory("fwlite_tt_lep_LO_order2_15weights_ref_delphes",  texName = "tt",      directory = [os.path.join( gen_dir, "fwlite_tt_lep_LO_order2_15weights_ref", "delphes" )], treeName = "Delphes") 
fwlite_WZ_lep_LO_order2_15weights_delphes           = Sample.fromDirectory("fwlite_WZ_lep_LO_order2_15weights_delphes",      texName = "WZ",      directory = [os.path.join( gen_dir, "fwlite_WZ_lep_LO_order2_15weights", "delphes" )], treeName = "Delphes")

#fwlite_tt_lep_LO_order2_15weights_ref.addFriend(  fwlite_tt_lep_LO_order2_15weights_ref_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_WZ_lep_LO_order2_15weights.addFriend(      fwlite_WZ_lep_LO_order2_15weights_delphes, treeName = "Delphes", sortFiles = True)


# hadronic + leptonic decays
fwlite_tt_LO_order2_15weights_ref               = Sample.fromDirectory("fwlite_tt_LO_order2_15weights_ref",  texName = "tt",      directory = [os.path.join( gen_dir, "fwlite_tt_LO_order2_15weights_ref" )]) 
fwlite_tt_LO_order2_15weights_ref.reweight_pkl  = gridpack_dir + "04062018/tt/order2/tt_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_tt_LO_order2_15weights_ref.nEvents       = 1000000 #? not checked!
fwlite_tt_LO_order2_15weights_ref.xsec          = 831.76 * 1610 / 499.9 #pb tt, xsec_SM_NNLO * xsec_BSM_LO / xsec_SM_LO

fwlite_WZ_LO_order2_15weights                   = Sample.fromDirectory("fwlite_WZ_LO_order2_15weights",     texName = "WZ",      directory = [os.path.join( gen_dir, "fwlite_WZ_LO_order2_15weights" )])
fwlite_WZ_LO_order2_15weights.reweight_pkl      = gridpack_dir + "04062018/WZ/order2/WZ_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_WZ_LO_order2_15weights.nEvents           = 1000000 #? not checked!
fwlite_WZ_LO_order2_15weights.xsec              = 47.13 #pb WZ, xsec_SM_NNLO

fwlite_tt_LO_order2_15weights_ref_delphes       = Sample.fromDirectory("fwlite_tt_LO_order2_15weights_ref_delphes",  texName = "tt",      directory = [os.path.join( gen_dir, "fwlite_tt_LO_order2_15weights_ref", "delphes" )], treeName = "Delphes") 
fwlite_WZ_LO_order2_15weights_delphes           = Sample.fromDirectory("fwlite_WZ_LO_order2_15weights_delphes",      texName = "WZ",      directory = [os.path.join( gen_dir, "fwlite_WZ_LO_order2_15weights", "delphes" )], treeName = "Delphes")

#fwlite_tt_LO_order2_15weights_ref.addFriend(  fwlite_tt_LO_order2_15weights_ref_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_WZ_LO_order2_15weights.addFriend(      fwlite_WZ_LO_order2_15weights_delphes, treeName = "Delphes", sortFiles = True)

"""
