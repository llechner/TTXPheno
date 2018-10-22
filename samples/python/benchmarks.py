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

##############################
##############################
##############################

# CMS

##############################
##############################
##############################

# signal samples with reference point

gen_dir = "/afs/hephy.at/data/llechner02/TTXPheno/skims/gen/v18/delphes_card_CMS/"

#fwlite_ttgammaLarge_LO_order2_15weights_ref_CMS              = Sample.fromDirectory("fwlite_ttgammaLarge_LO_order2_15weights_ref", texName = "ttgamma",  directory = [os.path.join( gen_dir, "fwlite_ttgammaLarge_LO_order2_15weights_ref" )])
#fwlite_ttgammaLarge_LO_order2_15weights_ref_CMS.reweight_pkl = gridpack_dir + "18052018_ref/ttgamma/order2/ttgamma0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
#fwlite_ttgammaLarge_LO_order2_15weights_ref_CMS.nEvents      = 10000000
#fwlite_ttgammaLarge_LO_order2_15weights_ref_CMS.xsec         = 7.092 * 3.697 / 2.302 #pb ttgamma gridpack * ttgamma NLO Daniel / ttgamma LO run.py UFO
#check

# reference point 15/2
fwlite_ttZ_ll_LO_order2_15weights_ref_CMS               = Sample.fromDirectory("fwlite_ttZ_ll_LO_order2_15weights_ref",  texName = "ttZ",      directory = [os.path.join( gen_dir, "fwlite_ttZ_ll_LO_order2_15weights_ref" )]) 
fwlite_ttZ_ll_LO_order2_15weights_ref_CMS.reweight_pkl  = gridpack_dir + "18052018_ref/ttZ/order2/ttZ0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttZ_ll_LO_order2_15weights_ref_CMS.nEvents       = 990000
fwlite_ttZ_ll_LO_order2_15weights_ref_CMS.xsec          = 0.5205 * 0.0915 / 0.0565 #pb ttZ, Z->ll, ttZ gridpack * ttZ NLO Daniel / ttZ LO run.py UFO
#check

fwlite_ttgamma_LO_order2_15weights_ref_CMS              = Sample.fromDirectory("fwlite_ttgamma_LO_order2_15weights_ref", texName = "ttgamma",  directory = [os.path.join( gen_dir, "fwlite_ttgamma_LO_order2_15weights_ref" )])
fwlite_ttgamma_LO_order2_15weights_ref_CMS.reweight_pkl = gridpack_dir + "18052018_ref/ttgamma/order2/ttgamma0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttgamma_LO_order2_15weights_ref_CMS.nEvents      = 970000
fwlite_ttgamma_LO_order2_15weights_ref_CMS.xsec         = 7.092 * 3.697 / 2.302 #pb ttgamma gridpack * ttgamma NLO Daniel / ttgamma LO run.py UFO
#check

#fwlite_ttgammaLarge_LO_order2_15weights_ref_CMS_delphes      = Sample.fromDirectory("fwlite_ttgammaLarge_LO_order2_15weights_ref_delphes", texName = "ttgamma",  directory = [os.path.join( gen_dir, "fwlite_ttgammaLarge_LO_order2_15weights_ref", "delphes" )], treeName = "Delphes")
#fwlite_ttZ_ll_LO_order2_15weights_ref_CMS_delphes       = Sample.fromDirectory("fwlite_ttZ_ll_LO_order2_15weights_ref_delphes",  texName = "ttZ",      directory = [os.path.join( gen_dir, "fwlite_ttZ_ll_LO_order2_15weights_ref", "delphes" )], treeName = "Delphes") 
#fwlite_ttW_LO_order2_15weights_ref_CMS_delphes          = Sample.fromDirectory("fwlite_ttW_LO_order2_15weights_ref_delphes",     texName = "ttW",      directory = [os.path.join( gen_dir, "fwlite_ttW_LO_order2_15weights_ref", "delphes" )], treeName = "Delphes")
#fwlite_ttgamma_LO_order2_15weights_ref_CMS_delphes      = Sample.fromDirectory("fwlite_ttgamma_LO_order2_15weights_ref_delphes", texName = "ttgamma",  directory = [os.path.join( gen_dir, "fwlite_ttgamma_LO_order2_15weights_ref", "delphes" )], treeName = "Delphes")

#fwlite_ttgammaLarge_LO_order2_15weights_ref_CMS.addFriend( fwlite_ttgammaLarge_LO_order2_15weights_ref_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_ttZ_ll_LO_order2_15weights_ref_CMS.addFriend(  fwlite_ttZ_ll_LO_order2_15weights_ref_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_ttW_LO_order2_15weights_ref_CMS.addFriend(     fwlite_ttW_LO_order2_15weights_ref_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_ttgamma_LO_order2_15weights_ref_CMS.addFriend( fwlite_ttgamma_LO_order2_15weights_ref_delphes, treeName = "Delphes", sortFiles = True)

##############################
##############################
##############################

# backgrounds without reference point
# leptonic decays W > lnu, Z > ll, t > Wb
fwlite_tt_dilep_LO_order2_15weights_CMS               = Sample.fromDirectory("fwlite_tt_dilep_LO_order2_15weights",  texName = "ttdilep",      directory = [os.path.join( gen_dir, "fwlite_tt_dilep_LO_order2_15weights" )]) 
fwlite_tt_dilep_LO_order2_15weights_CMS.reweight_pkl  = gridpack_dir + "06072018/tt/order2/tt_lep_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_tt_dilep_LO_order2_15weights_CMS.nEvents       = 1000000 
fwlite_tt_dilep_LO_order2_15weights_CMS.xsec          = 45.71 * 831.76 / 485.8 #tt dilep gridpack * tt NLO Daniel / tt LO run.py UFO
#check

fwlite_tt_nonhad_LO_order2_15weights_CMS               = Sample.fromDirectory("fwlite_tt_nonhad_LO_order2_15weights",  texName = "ttnonhad",      directory = [os.path.join( gen_dir, "fwlite_tt_nonhad_LO_order2_15weights" )]) 
fwlite_tt_nonhad_LO_order2_15weights_CMS.reweight_pkl  = gridpack_dir + "06072018/tt_semilep/order2/tt_semilep_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_tt_nonhad_LO_order2_15weights_CMS.nEvents       = 1000000 
fwlite_tt_nonhad_LO_order2_15weights_CMS.xsec          = 300.7 * 831.76 / 485.8 #tt nonhad gridpack * tt NLO Daniel / tt LO run.py UFO
#check

#fwlite_tt_full_LO_order2_15weights_CMS               = Sample.fromDirectory("fwlite_tt_full_LO_order2_15weights",  texName = "tt",      directory = [os.path.join( gen_dir, "fwlite_tt_full_LO_order2_15weights" )]) 
#fwlite_tt_full_LO_order2_15weights_CMS.reweight_pkl  = gridpack_dir + "06072018/tt_semilep/order2/tt_semilep_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
#fwlite_tt_full_LO_order2_15weights_CMS.nEvents       = 10000000 
#fwlite_tt_full_LO_order2_15weights_CMS.xsec          = 494.9 * 831.76 / 485.8 #tt full gridpack * tt NLO Daniel / tt LO run.py UFO
#check

#fwlite_tt_dilep_LO_order2_15weights_CMS_delphes     = Sample.fromDirectory("fwlite_tt_dilep_LO_order2_15weights_delphes",  texName = "ttdilep",      directory = [os.path.join( gen_dir, "fwlite_tt_dilep_LO_order2_15weights", "delphes" )], treeName = "Delphes") 
#fwlite_tt_nonhad_LO_order2_15weights_CMS_delphes    = Sample.fromDirectory("fwlite_tt_nonhad_LO_order2_15weights_delphes", texName = "ttnonhad",     directory = [os.path.join( gen_dir, "fwlite_tt_nonhad_LO_order2_15weights", "delphes" )], treeName = "Delphes") 
#fwlite_tt_full_LO_order2_15weights_CMS_delphes      = Sample.fromDirectory("fwlite_tt_full_LO_order2_15weights_delphes",   texName = "tt",           directory = [os.path.join( gen_dir, "fwlite_tt_full_LO_order2_15weights", "delphes" )], treeName = "Delphes") 

#fwlite_tt_full_LO_order2_15weights_CMS.addFriend(    fwlite_tt_full_LO_order2_15weights_CMS_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_tt_nonhad_LO_order2_15weights_CMS.addFriend(  fwlite_tt_nonhad_LO_order2_15weights_CMS_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_tt_dilep_LO_order2_15weights_CMS.addFriend(   fwlite_tt_dilep_LO_order2_15weights_CMS_delphes, treeName = "Delphes", sortFiles = True)

##############################
##############################
##############################

fwlite_tZq_LO_order2_15weights_CMS               = Sample.fromDirectory("fwlite_tZq_LO_order2_15weights",  texName = "tZq",      directory = [os.path.join( gen_dir, "fwlite_tZq_LO_order2_15weights" )]) 
fwlite_tZq_LO_order2_15weights_CMS.reweight_pkl  = gridpack_dir + "06072018/tZq/order2/tZq_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_tZq_LO_order2_15weights_CMS.nEvents       = 1000000 
fwlite_tZq_LO_order2_15weights_CMS.xsec          = 0.06774 * 0.0758 / 0.06148 #pb tZq gridpack * tZq NLO Daniel / tZq LO run.py DIM6
#check

fwlite_tW_LO_order2_15weights_CMS               = Sample.fromDirectory("fwlite_tW_LO_order2_15weights",  texName = "tW",      directory = [os.path.join( gen_dir, "fwlite_tW_LO_order2_15weights" )]) 
fwlite_tW_LO_order2_15weights_CMS.reweight_pkl  = gridpack_dir + "06072018/tW/order2/tW_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_tW_LO_order2_15weights_CMS.nEvents       = 1000000 
fwlite_tW_LO_order2_15weights_CMS.xsec          = 15.82 * 19.55 / 18.09 #pb tW gridpack * tW NLO Daniel / tW LO run.py UFO
#check

fwlite_tWZ_LO_order2_15weights_CMS               = Sample.fromDirectory("fwlite_tWZ_LO_order2_15weights",  texName = "tWZ",      directory = [os.path.join( gen_dir, "fwlite_tWZ_LO_order2_15weights" )]) 
fwlite_tWZ_LO_order2_15weights_CMS.reweight_pkl  = gridpack_dir + "06072018/tWZ/order2/tWZ_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_tWZ_LO_order2_15weights_CMS.nEvents       = 1000000 
fwlite_tWZ_LO_order2_15weights_CMS.xsec          = 0.01035 * 0.01123 / 0.01225 #tWZ gridpack * tWll NLO Daniel / tWll LO run.py UFO
#check

fwlite_Zgamma_LO_order2_15weights_CMS               = Sample.fromDirectory("fwlite_Zgamma_LO_order2_15weights",  texName = "Zgamma",      directory = [os.path.join( gen_dir, "fwlite_Zgamma_LO_order2_15weights" )]) 
fwlite_Zgamma_LO_order2_15weights_CMS.reweight_pkl  = gridpack_dir + "06072018/Zgamma/order2/Zgamma_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_Zgamma_LO_order2_15weights_CMS.nEvents       = 1000000 
fwlite_Zgamma_LO_order2_15weights_CMS.xsec          = 131.3 #pb ZGTo2LG_ext Daniel
#check?

fwlite_ttgamma_bg_LO_order2_15weights_CMS               = Sample.fromDirectory("fwlite_ttgamma_bg_LO_order2_15weights",  texName = "ttgamma",      directory = [os.path.join( gen_dir, "fwlite_ttgamma_bg_LO_order2_15weights" )]) 
fwlite_ttgamma_bg_LO_order2_15weights_CMS.reweight_pkl  = gridpack_dir + "06072018/ttgamma/order2/ttgamma_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttgamma_bg_LO_order2_15weights_CMS.nEvents       = 1000000 
fwlite_ttgamma_bg_LO_order2_15weights_CMS.xsec          = 2.179 * 3.697 / 2.302 #pb ttgamma gridpack * ttgamma NLO Daniel / ttgamma LO run.py UFO
#check

fwlite_WZ_lep_LO_order2_15weights_CMS                   = Sample.fromDirectory("fwlite_WZ_lep_LO_order2_15weights",     texName = "WZ",      directory = [os.path.join( gen_dir, "fwlite_WZ_lep_LO_order2_15weights" )])
fwlite_WZ_lep_LO_order2_15weights_CMS.reweight_pkl      = gridpack_dir + "06072018/WZ/order2/WZ_lep_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_WZ_lep_LO_order2_15weights_CMS.nEvents           = 1000000 
fwlite_WZ_lep_LO_order2_15weights_CMS.xsec              = 4.666 #WZTo3LNu_amcatnlo  #47.13*(3*0.108)*(3*0.0336) #pb WZ NLO Daniel * BR(Wlep) * BR(Zlep)
#check?

#fwlite_tZq_LO_order2_15weights_CMS_delphes           = Sample.fromDirectory("fwlite_tZq_LO_order2_15weights_delphes",         texName = "tZq",        directory = [os.path.join( gen_dir, "fwlite_tZq_LO_order2_15weights", "delphes" )], treeName = "Delphes") 
#fwlite_tW_LO_order2_15weights_CMS_delphes            = Sample.fromDirectory("fwlite_tW_LO_order2_15weights_delphes",          texName = "tW",         directory = [os.path.join( gen_dir, "fwlite_tW_LO_order2_15weights", "delphes" )], treeName = "Delphes") 
#fwlite_tWZ_LO_order2_15weights_CMS_delphes           = Sample.fromDirectory("fwlite_tWZ_LO_order2_15weights_delphes",         texName = "tWZ",        directory = [os.path.join( gen_dir, "fwlite_tWZ_LO_order2_15weights", "delphes" )], treeName = "Delphes") 
#fwlite_Zgamma_LO_order2_15weights_CMS_delphes        = Sample.fromDirectory("fwlite_Zgamma_LO_order2_15weights_delphes",      texName = "Zgamma",     directory = [os.path.join( gen_dir, "fwlite_Zgamma_LO_order2_15weights", "delphes" )], treeName = "Delphes") 
#fwlite_ttgamma_bg_LO_order2_15weights_CMS_delphes    = Sample.fromDirectory("fwlite_ttgamma_bg_LO_order2_15weights_delphes",  texName = "ttgamma_bg", directory = [os.path.join( gen_dir, "fwlite_ttgamma_bg_LO_order2_15weights", "delphes" )], treeName = "Delphes") 
#fwlite_WZ_lep_LO_order2_15weights_CMS_delphes        = Sample.fromDirectory("fwlite_WZ_lep_LO_order2_15weights_delphes",      texName = "WZ",      directory = [os.path.join( gen_dir, "fwlite_WZ_lep_LO_order2_15weights", "delphes" )], treeName = "Delphes")

#fwlite_tZq_LO_order2_15weights_CMS.addFriend(        fwlite_tZq_LO_order2_15weights_CMS_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_tW_LO_order2_15weights_CMS.addFriend(         fwlite_tW_LO_order2_15weights_CMS_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_tWZ_LO_order2_15weights_CMS.addFriend(        fwlite_tWz_LO_order2_15weights_CMS_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_Zgamma_LO_order2_15weights_CMS.addFriend(     fwlite_Zgamma_LO_order2_15weights_CMS_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_ttgamma_bg_LO_order2_15weights_CMS.addFriend( fwlite_ttgamma_bg_LO_order2_15weights_CMS_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_WZ_lep_LO_order2_15weights_CMS.addFriend(     fwlite_WZ_lep_LO_order2_15weights_CMS_delphes, treeName = "Delphes", sortFiles = True)

##############################
##############################
##############################

# CMS Phase2 200 PU

##############################
##############################
##############################

# signal samples with reference point

#gen_dir = "/afs/hephy.at/data/rschoefbeck01/TTXPheno/skims/gen/v22/"
gen_dir = "/afs/hephy.at/data/llechner02/TTXPheno/skims/gen/v23"

#fwlite_ttgammaLarge_LO_order2_15weights_ref_phase2_CMS              = Sample.fromDirectory("fwlite_ttgammaLarge_LO_order2_15weights_ref", texName = "ttgamma",  directory = [os.path.join( gen_dir, "fwlite_ttgammaLarge_LO_order2_15weights_ref" )])
#fwlite_ttgammaLarge_LO_order2_15weights_ref_phase2_CMS.reweight_pkl = gridpack_dir + "18052018_ref/ttgamma/order2/ttgamma0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
#fwlite_ttgammaLarge_LO_order2_15weights_ref_phase2_CMS.nEvents      = 10000000
#fwlite_ttgammaLarge_LO_order2_15weights_ref_phase2_CMS.xsec         = 7.092 * 3.697 / 2.302 #pb ttgamma gridpack * ttgamma NLO Daniel / ttgamma LO run.py UFO
#check

# reference point 15/2
fwlite_ttZ_ll_LO_order2_15weights_ref_phase2_CMS               = Sample.fromDirectory("fwlite_ttZ_ll_LO_order2_15weights_ref",  texName = "ttZ",      directory = [os.path.join( gen_dir, "fwlite_ttZ_ll_LO_order2_15weights_ref" )]) 
fwlite_ttZ_ll_LO_order2_15weights_ref_phase2_CMS.reweight_pkl  = gridpack_dir + "18052018_ref/ttZ/order2/ttZ0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttZ_ll_LO_order2_15weights_ref_phase2_CMS.nEvents       = 990000
fwlite_ttZ_ll_LO_order2_15weights_ref_phase2_CMS.xsec          = 0.5205 * (0.0915 / 0.0565) #pb ttZ, Z->ll, ttZ gridpack * ttZ NLO Daniel / ttZ LO run.py UFO
fwlite_ttZ_ll_LO_order2_15weights_ref_phase2_CMS.xsec14        = 0.5205 * (0.0915 / 0.0565) * (0.7152 / 0.616) #pb ttZ, Z->ll, ttZ gridpack * ttZ NLO Daniel / ttZ LO run.py UFO * ttZ jets 14 TeV / ttZ jets 13 TeV
#check
#16% 14 TeV correction

#fwlite_ttgamma_LO_order2_15weights_ref_phase2_CMS              = Sample.fromDirectory("fwlite_ttgamma_LO_order2_15weights_ref", texName = "ttgamma",  directory = [os.path.join( gen_dir, "fwlite_ttgamma_LO_order2_15weights_ref" )])
#fwlite_ttgamma_LO_order2_15weights_ref_phase2_CMS.reweight_pkl = gridpack_dir + "18052018_ref/ttgamma/order2/ttgamma0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
#fwlite_ttgamma_LO_order2_15weights_ref_phase2_CMS.nEvents      = 970000
#fwlite_ttgamma_LO_order2_15weights_ref_phase2_CMS.xsec         = 7.092 * 3.697 / 2.302 #pb ttgamma gridpack * ttgamma NLO Daniel / ttgamma LO run.py UFO
#check

#fwlite_ttgammaLarge_LO_order2_15weights_ref_phase2_CMS_delphes      = Sample.fromDirectory("fwlite_ttgammaLarge_LO_order2_15weights_ref_delphes", texName = "ttgamma",  directory = [os.path.join( gen_dir, "fwlite_ttgammaLarge_LO_order2_15weights_ref", "delphes" )], treeName = "Delphes")
#fwlite_ttZ_ll_LO_order2_15weights_ref_phase2_CMS_delphes       = Sample.fromDirectory("fwlite_ttZ_ll_LO_order2_15weights_ref_delphes",  texName = "ttZ",      directory = [os.path.join( gen_dir, "fwlite_ttZ_ll_LO_order2_15weights_ref", "delphes" )], treeName = "Delphes") 
#fwlite_ttW_LO_order2_15weights_ref_phase2_CMS_delphes          = Sample.fromDirectory("fwlite_ttW_LO_order2_15weights_ref_delphes",     texName = "ttW",      directory = [os.path.join( gen_dir, "fwlite_ttW_LO_order2_15weights_ref", "delphes" )], treeName = "Delphes")
#fwlite_ttgamma_LO_order2_15weights_ref_phase2_CMS_delphes      = Sample.fromDirectory("fwlite_ttgamma_LO_order2_15weights_ref_delphes", texName = "ttgamma",  directory = [os.path.join( gen_dir, "fwlite_ttgamma_LO_order2_15weights_ref", "delphes" )], treeName = "Delphes")

#fwlite_ttgammaLarge_LO_order2_15weights_ref_phase2_CMS.addFriend( fwlite_ttgammaLarge_LO_order2_15weights_ref_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_ttZ_ll_LO_order2_15weights_ref_phase2_CMS.addFriend(  fwlite_ttZ_ll_LO_order2_15weights_ref_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_ttW_LO_order2_15weights_ref_phase2_CMS.addFriend(     fwlite_ttW_LO_order2_15weights_ref_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_ttgamma_LO_order2_15weights_ref_phase2_CMS.addFriend( fwlite_ttgamma_LO_order2_15weights_ref_delphes, treeName = "Delphes", sortFiles = True)

##############################
##############################
##############################

# backgrounds without reference point
# leptonic decays W > lnu, Z > ll, t > Wb
#fwlite_tt_dilep_LO_order2_15weights_phase2_CMS               = Sample.fromDirectory("fwlite_tt_dilep_LO_order2_15weights",  texName = "ttdilep",      directory = [os.path.join( gen_dir, "fwlite_tt_dilep_LO_order2_15weights" )]) 
#fwlite_tt_dilep_LO_order2_15weights_phase2_CMS.reweight_pkl  = gridpack_dir + "06072018/tt/order2/tt_lep_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
#fwlite_tt_dilep_LO_order2_15weights_phase2_CMS.nEvents       = 1000000 
#fwlite_tt_dilep_LO_order2_15weights_phase2_CMS.xsec          = 45.71 * 831.76 / 485.8 #tt dilep gridpack * tt NLO Daniel / tt LO run.py UFO
#check

#fwlite_tt_nonhad_LO_order2_15weights_phase2_CMS               = Sample.fromDirectory("fwlite_tt_nonhad_LO_order2_15weights",  texName = "ttnonhad",      directory = [os.path.join( gen_dir, "fwlite_tt_nonhad_LO_order2_15weights" )]) 
#fwlite_tt_nonhad_LO_order2_15weights_phase2_CMS.reweight_pkl  = gridpack_dir + "06072018/tt_semilep/order2/tt_semilep_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
#fwlite_tt_nonhad_LO_order2_15weights_phase2_CMS.nEvents       = 1000000 
#fwlite_tt_nonhad_LO_order2_15weights_phase2_CMS.xsec          = 300.7 * 831.76 / 485.8 #tt nonhad gridpack * tt NLO Daniel / tt LO run.py UFO
#check

#fwlite_tt_full_LO_order2_15weights_phase2_CMS               = Sample.fromDirectory("fwlite_tt_full_LO_order2_15weights",  texName = "tt",      directory = [os.path.join( gen_dir, "fwlite_tt_full_LO_order2_15weights" )]) 
#fwlite_tt_full_LO_order2_15weights_phase2_CMS.reweight_pkl  = gridpack_dir + "06072018/tt_semilep/order2/tt_semilep_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
#fwlite_tt_full_LO_order2_15weights_phase2_CMS.nEvents       = 10000000 
#fwlite_tt_full_LO_order2_15weights_phase2_CMS.xsec          = 494.9 * 831.76 / 485.8 #tt full gridpack * tt NLO Daniel / tt LO run.py UFO
#check

#fwlite_tt_dilep_LO_order2_15weights_phase2_CMS_delphes     = Sample.fromDirectory("fwlite_tt_dilep_LO_order2_15weights_delphes",  texName = "ttdilep",      directory = [os.path.join( gen_dir, "fwlite_tt_dilep_LO_order2_15weights", "delphes" )], treeName = "Delphes") 
#fwlite_tt_nonhad_LO_order2_15weights_phase2_CMS_delphes    = Sample.fromDirectory("fwlite_tt_nonhad_LO_order2_15weights_delphes", texName = "ttnonhad",     directory = [os.path.join( gen_dir, "fwlite_tt_nonhad_LO_order2_15weights", "delphes" )], treeName = "Delphes") 
#fwlite_tt_full_LO_order2_15weights_phase2_CMS_delphes      = Sample.fromDirectory("fwlite_tt_full_LO_order2_15weights_delphes",   texName = "tt",           directory = [os.path.join( gen_dir, "fwlite_tt_full_LO_order2_15weights", "delphes" )], treeName = "Delphes") 

#fwlite_tt_full_LO_order2_15weights_phase2_CMS.addFriend(    fwlite_tt_full_LO_order2_15weights_phase2_CMS_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_tt_nonhad_LO_order2_15weights_phase2_CMS.addFriend(  fwlite_tt_nonhad_LO_order2_15weights_phase2_CMS_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_tt_dilep_LO_order2_15weights_phase2_CMS.addFriend(   fwlite_tt_dilep_LO_order2_15weights_phase2_CMS_delphes, treeName = "Delphes", sortFiles = True)

##############################
##############################
##############################

fwlite_tZq_LO_order2_15weights_phase2_CMS               = Sample.fromDirectory("fwlite_tZq_LO_order2_15weights",  texName = "tZq",      directory = [os.path.join( gen_dir, "fwlite_tZq_LO_order2_15weights" )]) 
fwlite_tZq_LO_order2_15weights_phase2_CMS.reweight_pkl  = gridpack_dir + "06072018/tZq/order2/tZq_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_tZq_LO_order2_15weights_phase2_CMS.nEvents       = 1000000 
fwlite_tZq_LO_order2_15weights_phase2_CMS.xsec          = 0.06774 * (0.0758 / 0.06148) #pb tZq gridpack * tZq NLO Daniel / tZq LO run.py DIM6
fwlite_tZq_LO_order2_15weights_phase2_CMS.xsec14        = 0.06774 * (0.0758 / 0.06148) * (0.085 / 0.0758) #pb tZq gridpack * tZq NLO Daniel / tZq LO run.py DIM6 * tZq ll 14TeV / tZq ll 13 TeV
#check
#12% 14 TeV correction

#fwlite_tW_LO_order2_15weights_phase2_CMS               = Sample.fromDirectory("fwlite_tW_LO_order2_15weights",  texName = "tW",      directory = [os.path.join( gen_dir, "fwlite_tW_LO_order2_15weights" )]) 
#fwlite_tW_LO_order2_15weights_phase2_CMS.reweight_pkl  = gridpack_dir + "06072018/tW/order2/tW_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
#fwlite_tW_LO_order2_15weights_phase2_CMS.nEvents       = 1000000 
#fwlite_tW_LO_order2_15weights_phase2_CMS.xsec          = 15.82 * 19.55 / 18.09 #pb tW gridpack * tW NLO Daniel / tW LO run.py UFO
#check

fwlite_tWZ_LO_order2_15weights_phase2_CMS               = Sample.fromDirectory("fwlite_tWZ_LO_order2_15weights",  texName = "tWZ",      directory = [os.path.join( gen_dir, "fwlite_tWZ_LO_order2_15weights" )]) 
fwlite_tWZ_LO_order2_15weights_phase2_CMS.reweight_pkl  = gridpack_dir + "06072018/tWZ/order2/tWZ_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_tWZ_LO_order2_15weights_phase2_CMS.nEvents       = 1000000 
fwlite_tWZ_LO_order2_15weights_phase2_CMS.xsec          = 0.01035 * (0.01123 / 0.01225) #tWZ gridpack * tWll NLO Daniel / tWll LO run.py UFO
fwlite_tWZ_LO_order2_15weights_phase2_CMS.xsec14        = 0.01035 * (0.01123 / 0.01225) * 1.12 #tWZ gridpack * tWll NLO Daniel / tWll LO run.py UFO #scaled with 12% similar to tZq, no SF found for 14TeV
#check

#fwlite_Zgamma_LO_order2_15weights_phase2_CMS               = Sample.fromDirectory("fwlite_Zgamma_LO_order2_15weights",  texName = "Zgamma",      directory = [os.path.join( gen_dir, "fwlite_Zgamma_LO_order2_15weights" )]) 
#fwlite_Zgamma_LO_order2_15weights_phase2_CMS.reweight_pkl  = gridpack_dir + "06072018/Zgamma/order2/Zgamma_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
#fwlite_Zgamma_LO_order2_15weights_phase2_CMS.nEvents       = 1000000 
#fwlite_Zgamma_LO_order2_15weights_phase2_CMS.xsec          = 131.3 #pb ZGTo2LG_ext Daniel
#check?

fwlite_ttgamma_bg_LO_order2_15weights_phase2_CMS               = Sample.fromDirectory("fwlite_ttgamma_bg_LO_order2_15weights",  texName = "ttgamma",      directory = [os.path.join( gen_dir, "fwlite_ttgamma_bg_LO_order2_15weights" )]) 
fwlite_ttgamma_bg_LO_order2_15weights_phase2_CMS.reweight_pkl  = gridpack_dir + "06072018/ttgamma/order2/ttgamma_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttgamma_bg_LO_order2_15weights_phase2_CMS.nEvents       = 1000000 
fwlite_ttgamma_bg_LO_order2_15weights_phase2_CMS.xsec          = 2.179 * (3.697 / 2.302) #pb ttgamma gridpack * ttgamma NLO Daniel / ttgamma LO run.py UFO
fwlite_ttgamma_bg_LO_order2_15weights_phase2_CMS.xsec14        = 2.179 * (3.697 / 2.302) * (0.6231 / 0.607)#pb ttgamma gridpack * ttgamma NLO Daniel / ttgamma LO run.py UFO
#check

fwlite_WZ_lep_LO_order2_15weights_phase2_CMS                   = Sample.fromDirectory("fwlite_WZ_lep_LO_order2_15weights",     texName = "WZ",      directory = [os.path.join( gen_dir, "fwlite_WZ_lep_LO_order2_15weights" )])
fwlite_WZ_lep_LO_order2_15weights_phase2_CMS.reweight_pkl      = gridpack_dir + "06072018/WZ/order2/WZ_lep_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_WZ_lep_LO_order2_15weights_phase2_CMS.nEvents           = 1000000 
fwlite_WZ_lep_LO_order2_15weights_phase2_CMS.xsec              = 4.666 #WZTo3LNu_amcatnlo  #47.13*(3*0.108)*(3*0.0336) #pb WZ NLO Daniel * BR(Wlep) * BR(Zlep)
fwlite_WZ_lep_LO_order2_15weights_phase2_CMS.xsec14             = 4.666 * (3.721 / 3.22) #WZTo3LNu_amcatnlo  #pb WZ NLO Daniel 13 TeV * ZZTo2l2Q 14 TeV / ZZTo2l2Q 13 TeV
#check?
#15.5% 14 TeV correction

#fwlite_tZq_LO_order2_15weights_phase2_CMS_delphes           = Sample.fromDirectory("fwlite_tZq_LO_order2_15weights_delphes",         texName = "tZq",        directory = [os.path.join( gen_dir, "fwlite_tZq_LO_order2_15weights", "delphes" )], treeName = "Delphes") 
#fwlite_tW_LO_order2_15weights_phase2_CMS_delphes            = Sample.fromDirectory("fwlite_tW_LO_order2_15weights_delphes",          texName = "tW",         directory = [os.path.join( gen_dir, "fwlite_tW_LO_order2_15weights", "delphes" )], treeName = "Delphes") 
#fwlite_tWZ_LO_order2_15weights_phase2_CMS_delphes           = Sample.fromDirectory("fwlite_tWZ_LO_order2_15weights_delphes",         texName = "tWZ",        directory = [os.path.join( gen_dir, "fwlite_tWZ_LO_order2_15weights", "delphes" )], treeName = "Delphes") 
#fwlite_Zgamma_LO_order2_15weights_phase2_CMS_delphes        = Sample.fromDirectory("fwlite_Zgamma_LO_order2_15weights_delphes",      texName = "Zgamma",     directory = [os.path.join( gen_dir, "fwlite_Zgamma_LO_order2_15weights", "delphes" )], treeName = "Delphes") 
#fwlite_ttgamma_bg_LO_order2_15weights_phase2_CMS_delphes    = Sample.fromDirectory("fwlite_ttgamma_bg_LO_order2_15weights_delphes",  texName = "ttgamma_bg", directory = [os.path.join( gen_dir, "fwlite_ttgamma_bg_LO_order2_15weights", "delphes" )], treeName = "Delphes") 
#fwlite_WZ_lep_LO_order2_15weights_phase2_CMS_delphes        = Sample.fromDirectory("fwlite_WZ_lep_LO_order2_15weights_delphes",      texName = "WZ",      directory = [os.path.join( gen_dir, "fwlite_WZ_lep_LO_order2_15weights", "delphes" )], treeName = "Delphes")

#fwlite_tZq_LO_order2_15weights_phase2_CMS.addFriend(        fwlite_tZq_LO_order2_15weights_phase2_CMS_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_tW_LO_order2_15weights_phase2_CMS.addFriend(         fwlite_tW_LO_order2_15weights_phase2_CMS_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_tWZ_LO_order2_15weights_phase2_CMS.addFriend(        fwlite_tWz_LO_order2_15weights_phase2_CMS_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_Zgamma_LO_order2_15weights_phase2_CMS.addFriend(     fwlite_Zgamma_LO_order2_15weights_phase2_CMS_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_ttgamma_bg_LO_order2_15weights_phase2_CMS.addFriend( fwlite_ttgamma_bg_LO_order2_15weights_phase2_CMS_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_WZ_lep_LO_order2_15weights_phase2_CMS.addFriend(     fwlite_WZ_lep_LO_order2_15weights_phase2_CMS_delphes, treeName = "Delphes", sortFiles = True)

