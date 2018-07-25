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

gridpack_dir = "/afs/hephy.at/data/llechner01/TTXPheno/gridpacks/"
gen_dir = "/afs/hephy.at/data/rschoefbeck02/TopEFT/skims/gen/v2/"

# Robert first ttZ_ll dim6top scan
dim6top_ttZ_ll_LO_highStat_scan                           = Sample.fromDirectory("fwlite_ttZ_ll_LO_highStat_scan", texName = "ttZ (scan)", directory = [os.path.join( gen_dir, "fwlite_ttZ_ll_LO_highStat_scan")])
dim6top_ttZ_ll_LO_highStat_scan.reweight_pkl              = "/afs/hephy.at/data/rschoefbeck02/TopEFT/results/gridpacks/ttZ0j_rwgt_patch_625_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"

dim6top_ttZ_ll_LO_currentplane_highStat_scan              = Sample.fromDirectory("dim6top_ttZ_ll_LO_currentplane_highStat_scan", texName = "ttZ (current scan)", directory = [os.path.join( gen_dir, "fwlite_ttZ_ll_LO_currentplane_highStat_scan")])
dim6top_ttZ_ll_LO_currentplane_highStat_scan.reweight_pkl = "/afs/hephy.at/data/rschoefbeck02/TopEFT/results/gridpacks/ttZ0j_rwgt_patch_currentplane_highStat_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"

test = Sample.fromFiles("test", files = ["/afs/hephy.at/data/rschoefbeck02/TTXPheno/skims/gen/v2/test/test.root"], texName = "test")
test.reweight_pkl = '/afs/cern.ch/user/l/llechner/public/gridpacks_data/order_3/ttZ/gridpacks/reweight_card.pkl'

'''
gen_dir = "/afs/hephy.at/data/rschoefbeck01/TTXPheno/skims/gen/v3/"

# no reference point 8/3
fwlite_ttZ_ll_LO_order3_8weights               = Sample.fromDirectory("fwlite_ttZ_ll_LO_order3_8weights",  texName = "ttZ",      directory = [os.path.join( gen_dir, "fwlite_ttZ_ll_LO_order3_8weights" )]) 
fwlite_ttZ_ll_LO_order3_8weights.reweight_pkl  = gridpack_dir + "07052018/ttZ/order3/ttZ0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttZ_ll_LO_order3_8weights.nEvents       = 995000
fwlite_ttZ_ll_LO_order3_8weights.xsec          = 0.0915 #pb ttZ, Z->ll
#xsec prob wrong

fwlite_ttW_LO_order3_8weights                  = Sample.fromDirectory("fwlite_ttW_LO_order3_8weights",     texName = "ttW",      directory = [os.path.join( gen_dir, "fwlite_ttW_LO_order3_8weights" )])
fwlite_ttW_LO_order3_8weights.reweight_pkl     = gridpack_dir + "07052018/ttW/order3/ttW0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttW_LO_order3_8weights.nEvents          = 1000000
fwlite_ttW_LO_order3_8weights.xsec             = 0.2043 #pb ttW, W->lnu
#xsec prob wrong

fwlite_ttgamma_LO_order3_8weights              = Sample.fromDirectory("fwlite_ttgamma_LO_order3_8weights", texName = "ttgamma",  directory = [os.path.join( gen_dir, "fwlite_ttgamma_LO_order3_8weights" )])
fwlite_ttgamma_LO_order3_8weights.reweight_pkl = gridpack_dir + "07052018/ttgamma/order3/ttgamma0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttgamma_LO_order3_8weights.nEvents      = 1000000
fwlite_ttgamma_LO_order3_8weights.xsec         = 3.697 #pb ttgamma
#xsec prob wrong

fwlite_ttZ_ll_LO_order3_8weights_delphes       = Sample.fromDirectory("fwlite_ttZ_ll_LO_order3_8weights_delphes",  texName = "ttZ",      directory = [os.path.join( gen_dir, "fwlite_ttZ_ll_LO_order3_8weights", "delphes" )], treeName = "Delphes") 
fwlite_ttW_LO_order3_8weights_delphes          = Sample.fromDirectory("fwlite_ttW_LO_order3_8weights_delphes",     texName = "ttW",      directory = [os.path.join( gen_dir, "fwlite_ttW_LO_order3_8weights", "delphes" )], treeName = "Delphes")
fwlite_ttgamma_LO_order3_8weights_delphes      = Sample.fromDirectory("fwlite_ttgamma_LO_order3_8weights_delphes", texName = "ttgamma",  directory = [os.path.join( gen_dir, "fwlite_ttgamma_LO_order3_8weights", "delphes" )], treeName = "Delphes")

#fwlite_ttZ_ll_LO_order3_8weights.addFriend(  fwlite_ttZ_ll_LO_order3_8weights_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_ttW_LO_order3_8weights.addFriend(     fwlite_ttW_LO_order3_8weights_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_ttgamma_LO_order3_8weights.addFriend( fwlite_ttgamma_LO_order3_8weights_delphes, treeName = "Delphes", sortFiles = True)

gen_dir = "/afs/hephy.at/data/llechner01/TTXPheno/skims/gen/v3/"

# no reference point 15/2
fwlite_ttZ_ll_LO_order2_15weights               = Sample.fromDirectory("fwlite_ttZ_ll_LO_order2_15weights",  texName = "ttZ",      directory = [os.path.join( gen_dir, "fwlite_ttZ_ll_LO_order2_15weights" )]) 
fwlite_ttZ_ll_LO_order2_15weights.reweight_pkl  = gridpack_dir + "18052018/ttZ/order2/ttZ0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttZ_ll_LO_order2_15weights.nEvents       = 975000
fwlite_ttZ_ll_LO_order2_15weights.xsec          = 0.0915 #pb ttZ, Z->ll
#xsec prob wrong

fwlite_ttW_LO_order2_15weights                  = Sample.fromDirectory("fwlite_ttW_LO_order2_15weights",     texName = "ttW",      directory = [os.path.join( gen_dir, "fwlite_ttW_LO_order2_15weights" )])
fwlite_ttW_LO_order2_15weights.reweight_pkl     = gridpack_dir + "18052018/ttW/order2/ttW0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttW_LO_order2_15weights.nEvents          = 1000000
fwlite_ttW_LO_order2_15weights.xsec             = 0.2043 #pb ttW, W->lnu
#xsec prob wrong

fwlite_ttgamma_LO_order2_15weights              = Sample.fromDirectory("fwlite_ttgamma_LO_order2_15weights", texName = "ttgamma",  directory = [os.path.join( gen_dir, "fwlite_ttgamma_LO_order2_15weights" )])
fwlite_ttgamma_LO_order2_15weights.reweight_pkl = gridpack_dir + "18052018/ttgamma/order2/ttgamma0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttgamma_LO_order2_15weights.nEvents      = 990000
fwlite_ttgamma_LO_order2_15weights.xsec         = 3.697 #pb ttgamma
#xsec prob wrong

fwlite_ttZ_ll_LO_order2_15weights_delphes       = Sample.fromDirectory("fwlite_ttZ_ll_LO_order2_15weights_delphes",  texName = "ttZ",      directory = [os.path.join( gen_dir, "fwlite_ttZ_ll_LO_order2_15weights", "delphes" )], treeName = "Delphes") 
fwlite_ttW_LO_order2_15weights_delphes          = Sample.fromDirectory("fwlite_ttW_LO_order2_15weights_delphes",     texName = "ttW",      directory = [os.path.join( gen_dir, "fwlite_ttW_LO_order2_15weights", "delphes" )], treeName = "Delphes")
fwlite_ttgamma_LO_order2_15weights_delphes      = Sample.fromDirectory("fwlite_ttgamma_LO_order2_15weights_delphes", texName = "ttgamma",  directory = [os.path.join( gen_dir, "fwlite_ttgamma_LO_order2_15weights", "delphes" )], treeName = "Delphes")

#fwlite_ttZ_ll_LO_order2_15weights.addFriend(  fwlite_ttZ_ll_LO_order2_15weights_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_ttW_LO_order2_15weights.addFriend(     fwlite_ttW_LO_order2_15weights_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_ttgamma_LO_order2_15weights.addFriend( fwlite_ttgamma_LO_order2_15weights_delphes, treeName = "Delphes", sortFiles = True)

fwlite_ttZ_ll_LO_order2_15weights_ref_old               = Sample.fromDirectory("fwlite_ttZ_ll_LO_order2_15weights_ref",  texName = "ttZ",      directory = [os.path.join( gen_dir, "fwlite_ttZ_ll_LO_order2_15weights_ref" )]) 
fwlite_ttZ_ll_LO_order2_15weights_ref_old.reweight_pkl  = gridpack_dir + "18052018_ref/ttZ/order2/ttZ0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttZ_ll_LO_order2_15weights_ref_old.nEvents       = 990000
fwlite_ttZ_ll_LO_order2_15weights_ref_old.xsec          = 0.0915 * 0.6823 / 0.07082 #pb ttZ, Z->ll, xsec_SM_NNLO * xsec_BSM_LO / xsec_SM_LO
#xsec prob wrong

fwlite_ttgamma_LO_order2_15weights_ref_old              = Sample.fromDirectory("fwlite_ttgamma_LO_order2_15weights_ref", texName = "ttgamma",  directory = [os.path.join( gen_dir, "fwlite_ttgamma_LO_order2_15weights_ref" )])
fwlite_ttgamma_LO_order2_15weights_ref_old.reweight_pkl = gridpack_dir + "18052018_ref/ttgamma/order2/ttgamma0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttgamma_LO_order2_15weights_ref_old.nEvents      = 970000
fwlite_ttgamma_LO_order2_15weights_ref_old.xsec         = 3.697 * 7.838 / 2.439 #pb ttgamma, xsec_SM_NNLO * xsec_BSM_LO / xsec_SM_LO
#xsec prob wrong

fwlite_ttW_LO_order2_15weights_ref_old                  = Sample.fromDirectory("fwlite_ttW_LO_order2_15weights_ref",     texName = "ttW",      directory = [os.path.join( gen_dir, "fwlite_ttW_LO_order2_15weights_ref" )])
fwlite_ttW_LO_order2_15weights_ref_old.reweight_pkl     = gridpack_dir + "18052018_ref/ttW/order2/ttW0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttW_LO_order2_15weights_ref_old.nEvents          = 945000
fwlite_ttW_LO_order2_15weights_ref_old.xsec             = 0.2043 * 0.4097 / 0.1323 #pb ttW, W->lnu, xsec_SM_NNLO * xsec_BSM_LO / xsec_SM_LO
#xsec prob wrong

'''
gen_dir = "/afs/hephy.at/data/rschoefbeck01/TTXPheno/skims/gen/v13/"
#gen_dir = "/afs/hephy.at/data/llechner01/TTXPheno/skims/gen/v12/"

# reference point 15/2
fwlite_ttZ_ll_LO_order2_15weights_ref               = Sample.fromDirectory("fwlite_ttZ_ll_LO_order2_15weights_ref",  texName = "ttZ",      directory = [os.path.join( gen_dir, "fwlite_ttZ_ll_LO_order2_15weights_ref" )]) 
fwlite_ttZ_ll_LO_order2_15weights_ref.reweight_pkl  = gridpack_dir + "18052018_ref/ttZ/order2/ttZ0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttZ_ll_LO_order2_15weights_ref.nEvents       = 990000
fwlite_ttZ_ll_LO_order2_15weights_ref.xsec          = 0.5205 * 0.0915 / 0.07082 #pb ttZ, LO BSM (gridpack) * NLO SM / LO SM

fwlite_ttW_LO_order2_15weights_ref                  = Sample.fromDirectory("fwlite_ttW_LO_order2_15weights_ref",     texName = "ttW",      directory = [os.path.join( gen_dir, "fwlite_ttW_LO_order2_15weights_ref" )])
fwlite_ttW_LO_order2_15weights_ref.reweight_pkl     = gridpack_dir + "18052018_ref/ttW/order2/ttW0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttW_LO_order2_15weights_ref.nEvents          = 945000
fwlite_ttW_LO_order2_15weights_ref.xsec             = 0.3599 * 0.2043 / 0.1323 #pb ttW, LO BSM (gridpack) * NLO SM / LO SM

fwlite_ttgamma_LO_order2_15weights_ref              = Sample.fromDirectory("fwlite_ttgamma_LO_order2_15weights_ref", texName = "ttgamma",  directory = [os.path.join( gen_dir, "fwlite_ttgamma_LO_order2_15weights_ref" )])
fwlite_ttgamma_LO_order2_15weights_ref.reweight_pkl = gridpack_dir + "18052018_ref/ttgamma/order2/ttgamma0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttgamma_LO_order2_15weights_ref.nEvents      = 970000
fwlite_ttgamma_LO_order2_15weights_ref.xsec         = 7.092 * 3.697 / 2.439 #pb ttgamma, LO BSM (gridpack) * NLO SM / LO SM
#wrong sample

"""
fwlite_ttgamma1l_LO_order2_15weights_ref              = Sample.fromDirectory("fwlite_ttgamma1l_LO_order2_15weights_ref", texName = "ttgamma1l",  directory = [os.path.join( gen_dir, "fwlite_ttgamma1l_LO_order2_15weights_ref" )])
fwlite_ttgamma1l_LO_order2_15weights_ref.reweight_pkl = gridpack_dir + "24072018/ttgamma1l/order2/ttgamma1l_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttgamma1l_LO_order2_15weights_ref.nEvents      = 1000000
#fwlite_ttgamma1l_LO_order2_15weights_ref.xsec         = * 3.697 / 2.439 #pb ttgamma, xsec_SM_NNLO * xsec_BSM_LO / xsec_SM_LO
#check?

fwlite_ttgamma2l_LO_order2_15weights_ref              = Sample.fromDirectory("fwlite_ttgamma2l_LO_order2_15weights_ref", texName = "ttgamma2l",  directory = [os.path.join( gen_dir, "fwlite_ttgamma2l_LO_order2_15weights_ref" )])
fwlite_ttgamma2l_LO_order2_15weights_ref.reweight_pkl = gridpack_dir + "24072018/ttgamma2l/order2/ttgamma2l_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttgamma2l_LO_order2_15weights_ref.nEvents      = 1000000
#fwlite_ttgamma2l_LO_order2_15weights_ref.xsec         =  * 3.697 / 2.439 #pb ttgamma, xsec_SM_NNLO * xsec_BSM_LO / xsec_SM_LO
#check?
"""

fwlite_ttZ_ll_LO_order2_15weights_ref_delphes       = Sample.fromDirectory("fwlite_ttZ_ll_LO_order2_15weights_ref_delphes",  texName = "ttZ",      directory = [os.path.join( gen_dir, "fwlite_ttZ_ll_LO_order2_15weights_ref", "delphes" )], treeName = "Delphes") 
fwlite_ttW_LO_order2_15weights_ref_delphes          = Sample.fromDirectory("fwlite_ttW_LO_order2_15weights_ref_delphes",     texName = "ttW",      directory = [os.path.join( gen_dir, "fwlite_ttW_LO_order2_15weights_ref", "delphes" )], treeName = "Delphes")
fwlite_ttgamma_LO_order2_15weights_ref_delphes      = Sample.fromDirectory("fwlite_ttgamma_LO_order2_15weights_ref_delphes", texName = "ttgamma",  directory = [os.path.join( gen_dir, "fwlite_ttgamma_LO_order2_15weights_ref", "delphes" )], treeName = "Delphes")
#fwlite_ttgamma1l_LO_order2_15weights_ref_delphes      = Sample.fromDirectory("fwlite_ttgamma1l_LO_order2_15weights_ref_delphes", texName = "ttgamma1l",  directory = [os.path.join( gen_dir, "fwlite_ttgamma1l_LO_order2_15weights_ref", "delphes" )], treeName = "Delphes")
#fwlite_ttgamma2l_LO_order2_15weights_ref_delphes      = Sample.fromDirectory("fwlite_ttgamma2l_LO_order2_15weights_ref_delphes", texName = "ttgamma2l",  directory = [os.path.join( gen_dir, "fwlite_ttgamma2l_LO_order2_15weights_ref", "delphes" )], treeName = "Delphes")

#fwlite_ttZ_ll_LO_order2_15weights_ref.addFriend(  fwlite_ttZ_ll_LO_order2_15weights_ref_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_ttW_LO_order2_15weights_ref.addFriend(     fwlite_ttW_LO_order2_15weights_ref_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_ttgamma_LO_order2_15weights_ref.addFriend( fwlite_ttgamma_LO_order2_15weights_ref_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_ttgamma1l_LO_order2_15weights_ref.addFriend( fwlite_ttgamma1l_LO_order2_15weights_ref_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_ttgamma2l_LO_order2_15weights_ref.addFriend( fwlite_ttgamma2l_LO_order2_15weights_ref_delphes, treeName = "Delphes", sortFiles = True)

# backgrounds 15/2 (tt) or 0/2 (WZ)
# leptonic decays W > lnu, Z > ll, t > Wb
fwlite_tt_lep_LO_order2_15weights               = Sample.fromDirectory("fwlite_tt_lep_LO_order2_15weights",  texName = "tt",      directory = [os.path.join( gen_dir, "fwlite_tt_lep_LO_order2_15weights" )]) 
fwlite_tt_lep_LO_order2_15weights.reweight_pkl  = gridpack_dir + "06072018/tt/order2/tt_lep_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_tt_lep_LO_order2_15weights.nEvents       = 1000000 #? not checked!
fwlite_tt_lep_LO_order2_15weights.xsec          = 45.71*831.76/499.9 #pb LO (gridpack) * TTLep_pow NLO / LO, tt, proc card: p p > t t~, (t > w+ b, w+ > l+ vl), (t~ > w- b~, w- > l- vl~)

fwlite_tt_semilep_LO_order2_15weights               = Sample.fromDirectory("fwlite_tt_semilep_LO_order2_15weights",  texName = "ttsemilep",      directory = [os.path.join( gen_dir, "fwlite_tt_semilep_LO_order2_15weights" )]) 
fwlite_tt_semilep_LO_order2_15weights.reweight_pkl  = gridpack_dir + "06072018/tt_semilep/order2/tt_semilep_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_tt_semilep_LO_order2_15weights.nEvents       = 1000000 #? not checked!
fwlite_tt_semilep_LO_order2_15weights.xsec          = 300.7*831.76/499.9 #pb LO (gridpack) * TTSemiLep_pow / LO , tt, proc card: p p > t t~, (t~ > w- b~, w- > l- vl~) +  p p > t t~, (t > w+ b, w+ > l+ vl)

"""
fwlite_tt1l_LO_order2_15weights               = Sample.fromDirectory("fwlite_tt1l_LO_order2_15weights",  texName = "tt1l",      directory = [os.path.join( gen_dir, "fwlite_tt1l_LO_order2_15weights" )]) 
fwlite_tt1l_LO_order2_15weights.reweight_pkl  = gridpack_dir + "24072018/tt1l/order2/tt1l_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_tt1l_LO_order2_15weights.nEvents       = 1000000 #? not checked!
fwlite_tt1l_LO_order2_15weights.xsec          =  * 831.76 / 499.9 #pb tt1l*tt_NNLO/tt_LO, tt, proc card: p p > t t~, (t > w+ b, w+ > l+ vl), (t~ > w- b~, w- > l- vl~)

fwlite_tt2l_LO_order2_15weights               = Sample.fromDirectory("fwlite_tt2l_LO_order2_15weights",  texName = "tt2l",      directory = [os.path.join( gen_dir, "fwlite_tt2l_LO_order2_15weights" )]) 
fwlite_tt2l_LO_order2_15weights.reweight_pkl  = gridpack_dir + "24072018/tt2l/order2/tt2l_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_tt2l_LO_order2_15weights.nEvents       = 1000000 #? not checked!
fwlite_tt2l_LO_order2_15weights.xsec          =  * 831.76 / 499.9 #pb tt2l*tt_NNLO/tt_LO, tt, proc card: p p > t t~, (t~ > w- b~, w- > l- vl~) +  p p > t t~, (t > w+ b, w+ > l+ vl)
"""

fwlite_tZq_LO_order2_15weights               = Sample.fromDirectory("fwlite_tZq_LO_order2_15weights",  texName = "tZq",      directory = [os.path.join( gen_dir, "fwlite_tZq_LO_order2_15weights" )]) 
fwlite_tZq_LO_order2_15weights.reweight_pkl  = gridpack_dir + "06072018/tZq/order2/tZq_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_tZq_LO_order2_15weights.nEvents       = 1000000 #? not checked!
fwlite_tZq_LO_order2_15weights.xsec          = 0.06774 * 0.09418 / 0.06148 #pb LO * tZq_ll_NNLO / tZq_ll_LO, proc card: p p > t j l+ l-,  p p > t~ j l+ l-

fwlite_tW_LO_order2_15weights               = Sample.fromDirectory("fwlite_tW_LO_order2_15weights",  texName = "tW",      directory = [os.path.join( gen_dir, "fwlite_tW_LO_order2_15weights" )]) 
fwlite_tW_LO_order2_15weights.reweight_pkl  = gridpack_dir + "06072018/tW/order2/tW_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_tW_LO_order2_15weights.nEvents       = 1000000 #? not checked!
fwlite_tW_LO_order2_15weights.xsec          = 15.82 * 19.55 / 17.62 #pb LO BSM (gridpack) * NLO SM / LO SM, T_tWch_noFullyHad (x2 for t and tbar, /2 as here ether W or t decays leptonically (double counting)), proc card: p p > t l- vl~, p p > t~ l+ vl

fwlite_tWZ_LO_order2_15weights               = Sample.fromDirectory("fwlite_tWZ_LO_order2_15weights",  texName = "tWZ",      directory = [os.path.join( gen_dir, "fwlite_tWZ_LO_order2_15weights" )]) 
fwlite_tWZ_LO_order2_15weights.reweight_pkl  = gridpack_dir + "06072018/tWZ/order2/tWZ_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_tWZ_LO_order2_15weights.nEvents       = 1000000 #? not checked!
fwlite_tWZ_LO_order2_15weights.xsec          = 0.01035 * 0.01123 / 0.01191 #pb LO (gridpack) * tWll_NNLO / tWll_LO, proc card: p p > t w l+ l-, p p > t~ w l+ l-, wdecay = j l+ vl l- vl~

fwlite_Zgamma_LO_order2_15weights               = Sample.fromDirectory("fwlite_Zgamma_LO_order2_15weights",  texName = "Zgamma",      directory = [os.path.join( gen_dir, "fwlite_Zgamma_LO_order2_15weights" )]) 
fwlite_Zgamma_LO_order2_15weights.reweight_pkl  = gridpack_dir + "06072018/Zgamma/order2/Zgamma_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_Zgamma_LO_order2_15weights.nEvents       = 1000000 #? not checked!
fwlite_Zgamma_LO_order2_15weights.xsec          = 131.3 #pb ZGTo2LG_ext, proc card: p p > l+ l- a + p p > l+ l- a j
#check?

fwlite_ttgamma_bg_LO_order2_15weights               = Sample.fromDirectory("fwlite_ttgamma_bg_LO_order2_15weights",  texName = "ttgamma",      directory = [os.path.join( gen_dir, "fwlite_ttgamma_bg_LO_order2_15weights" )]) 
fwlite_ttgamma_bg_LO_order2_15weights.reweight_pkl  = gridpack_dir + "06072018/ttgamma/order2/ttgamma_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttgamma_bg_LO_order2_15weights.nEvents       = 1000000 #? not checked!
fwlite_ttgamma_bg_LO_order2_15weights.xsec          = 2.179 * 3.697 / 2.439 #pb LO (gridpack) * NLO TTGJets_ext / LO TTGJets_ext, ttgamma, proc card: p p > t t~ a

fwlite_WZ_lep_LO_order2_15weights                   = Sample.fromDirectory("fwlite_WZ_lep_LO_order2_15weights",     texName = "WZ",      directory = [os.path.join( gen_dir, "fwlite_WZ_lep_LO_order2_15weights" )])
fwlite_WZ_lep_LO_order2_15weights.reweight_pkl      = gridpack_dir + "06072018/WZ/order2/WZ_lep_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_WZ_lep_LO_order2_15weights.nEvents           = 1000000 #? not checked!
fwlite_WZ_lep_LO_order2_15weights.xsec              = 47.13*(3*0.108)*2*(3*0.0336) #pb WZ*BR(Wlep)*2(W+W-)*BR(Zlep), proc card: p p > l+ vl l- l+, p p > l- vl~ l- l+
#check?

fwlite_tt_lep_LO_order2_15weights_delphes        = Sample.fromDirectory("fwlite_tt_lep_LO_order2_15weights_ref_delphes",  texName = "tt",      directory = [os.path.join( gen_dir, "fwlite_tt_lep_LO_order2_15weights", "delphes" )], treeName = "Delphes") 
#fwlite_tt1l_LO_order2_15weights_delphes        = Sample.fromDirectory("fwlite_tt1l_LO_order2_15weights_ref_delphes",  texName = "tt1l",      directory = [os.path.join( gen_dir, "fwlite_tt1l_LO_order2_15weights", "delphes" )], treeName = "Delphes") 
#fwlite_tt2l_LO_order2_15weights_delphes        = Sample.fromDirectory("fwlite_tt2l_LO_order2_15weights_ref_delphes",  texName = "tt2l",      directory = [os.path.join( gen_dir, "fwlite_tt2l_LO_order2_15weights", "delphes" )], treeName = "Delphes") 
fwlite_tt_semilep_LO_order2_15weights_delphes    = Sample.fromDirectory("fwlite_tt_semilep_LO_order2_15weights_ref_delphes",  texName = "tt",         directory = [os.path.join( gen_dir, "fwlite_tt_semilep_LO_order2_15weights", "delphes" )], treeName = "Delphes") 
fwlite_tZq_LO_order2_15weights_delphes           = Sample.fromDirectory("fwlite_tZq_LO_order2_15weights_ref_delphes",         texName = "tZq",        directory = [os.path.join( gen_dir, "fwlite_tZq_LO_order2_15weights", "delphes" )], treeName = "Delphes") 
fwlite_tW_LO_order2_15weights_delphes            = Sample.fromDirectory("fwlite_tW_LO_order2_15weights_ref_delphes",          texName = "tW",         directory = [os.path.join( gen_dir, "fwlite_tW_LO_order2_15weights", "delphes" )], treeName = "Delphes") 
fwlite_tWZ_LO_order2_15weights_delphes           = Sample.fromDirectory("fwlite_tWZ_LO_order2_15weights_ref_delphes",         texName = "tWZ",        directory = [os.path.join( gen_dir, "fwlite_tWZ_LO_order2_15weights", "delphes" )], treeName = "Delphes") 
fwlite_Zgamma_LO_order2_15weights_delphes        = Sample.fromDirectory("fwlite_Zgamma_LO_order2_15weights_ref_delphes",      texName = "Zgamma",     directory = [os.path.join( gen_dir, "fwlite_Zgamma_LO_order2_15weights", "delphes" )], treeName = "Delphes") 
fwlite_ttgamma_bg_LO_order2_15weights_delphes    = Sample.fromDirectory("fwlite_ttgamma_bg_LO_order2_15weights_ref_delphes",  texName = "ttgamma_bg", directory = [os.path.join( gen_dir, "fwlite_ttgamma_bg_LO_order2_15weights", "delphes" )], treeName = "Delphes") 
fwlite_WZ_lep_LO_order2_15weights_delphes        = Sample.fromDirectory("fwlite_WZ_lep_LO_order2_15weights_delphes",      texName = "WZ",      directory = [os.path.join( gen_dir, "fwlite_WZ_lep_LO_order2_15weights", "delphes" )], treeName = "Delphes")

#fwlite_tt_lep_LO_order2_15weights.addFriend(     fwlite_tt_lep_LO_order2_15weights_ref_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_tt1l_LO_order2_15weights.addFriend(       fwlite_tt1l_LO_order2_15weights_ref_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_tt2l_LO_order2_15weights.addFriend(       fwlite_tt2l_LO_order2_15weights_ref_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_tt_semilep_LO_order2_15weights.addFriend( fwlite_tt_semilep_LO_order2_15weights_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_tZq_LO_order2_15weights.addFriend(        fwlite_tZq_LO_order2_15weights_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_tW_LO_order2_15weights.addFriend(         fwlite_tW_LO_order2_15weights_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_tWZ_LO_order2_15weights.addFriend(        fwlite_tWz_LO_order2_15weights_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_Zgamma_LO_order2_15weights.addFriend(     fwlite_Zgamma_LO_order2_15weights_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_ttgamma_bg_LO_order2_15weights.addFriend( fwlite_ttgamma_bg_LO_order2_15weights_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_WZ_lep_LO_order2_15weights.addFriend(     fwlite_WZ_lep_LO_order2_15weights_delphes, treeName = "Delphes", sortFiles = True)



"""
# backgrounds 15/2 (tt) or 0/2 (WZ)
# leptonic decays W > lnu, Z > ll, t > Wb
fwlite_tt_lep_LO_order2_15weights_ref               = Sample.fromDirectory("fwlite_tt_lep_LO_order2_15weights_ref",  texName = "tt",      directory = [os.path.join( gen_dir, "fwlite_tt_lep_LO_order2_15weights_ref" )]) 
fwlite_tt_lep_LO_order2_15weights_ref.reweight_pkl  = gridpack_dir + "05062018_ref/tt/order2/tt_lep_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_tt_lep_LO_order2_15weights_ref.nEvents       = 1000000 #? not checked!
fwlite_tt_lep_LO_order2_15weights_ref.xsec          = 87.3148 * 735.1 / 46.57 #pb tt, W->lnu, xsec_SM_NNLO * xsec_BSM_LO / xsec_SM_LO
#xsec prob wrong

# hadronic + leptonic decays
fwlite_tt_LO_order2_15weights_ref               = Sample.fromDirectory("fwlite_tt_LO_order2_15weights_ref",  texName = "tt",      directory = [os.path.join( gen_dir, "fwlite_tt_LO_order2_15weights_ref" )]) 
fwlite_tt_LO_order2_15weights_ref.reweight_pkl  = gridpack_dir + "04062018/tt/order2/tt_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_tt_LO_order2_15weights_ref.nEvents       = 1000000 #? not checked!
fwlite_tt_LO_order2_15weights_ref.xsec          = 831.76 * 1610 / 499.9 #pb tt, xsec_SM_NNLO * xsec_BSM_LO / xsec_SM_LO
#xsec prob wrong

fwlite_WZ_LO_order2_15weights                   = Sample.fromDirectory("fwlite_WZ_LO_order2_15weights",     texName = "WZ",      directory = [os.path.join( gen_dir, "fwlite_WZ_LO_order2_15weights" )])
fwlite_WZ_LO_order2_15weights.reweight_pkl      = gridpack_dir + "04062018/WZ/order2/WZ_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_WZ_LO_order2_15weights.nEvents           = 1000000 #? not checked!
fwlite_WZ_LO_order2_15weights.xsec              = 47.13 #pb WZ, xsec_SM_NNLO
#xsec prob wrong

fwlite_tt_lep_LO_order2_15weights_ref_delphes   = Sample.fromDirectory("fwlite_tt_lep_LO_order2_15weights_ref_delphes",  texName = "tt",      directory = [os.path.join( gen_dir, "fwlite_tt_lep_LO_order2_15weights_ref", "delphes" )], treeName = "Delphes") 
fwlite_tt_LO_order2_15weights_ref_delphes       = Sample.fromDirectory("fwlite_tt_LO_order2_15weights_ref_delphes",  texName = "tt",      directory = [os.path.join( gen_dir, "fwlite_tt_LO_order2_15weights_ref", "delphes" )], treeName = "Delphes") 
fwlite_WZ_LO_order2_15weights_delphes           = Sample.fromDirectory("fwlite_WZ_LO_order2_15weights_delphes",      texName = "WZ",      directory = [os.path.join( gen_dir, "fwlite_WZ_LO_order2_15weights", "delphes" )], treeName = "Delphes")

#fwlite_tt_lep_LO_order2_15weights_ref.addFriend(  fwlite_tt_lep_LO_order2_15weights_ref_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_tt_LO_order2_15weights_ref.addFriend(  fwlite_tt_LO_order2_15weights_ref_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_WZ_LO_order2_15weights.addFriend(      fwlite_WZ_LO_order2_15weights_delphes, treeName = "Delphes", sortFiles = True)
"""
