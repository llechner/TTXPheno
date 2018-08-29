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

# signal samples with reference point

#gen_dir = "/afs/hephy.at/data/llechner02/TTXPheno/skims/gen/v17/"
gen_dir = "/afs/hephy.at/data/rschoefbeck01/TTXPheno/skims/gen/v18/"

fwlite_ttgammaLarge_LO_order2_15weights_ref              = Sample.fromDirectory("fwlite_ttgammaLarge_LO_order2_15weights_ref", texName = "ttgamma",  directory = [os.path.join( gen_dir, "fwlite_ttgammaLarge_LO_order2_15weights_ref" )])
fwlite_ttgammaLarge_LO_order2_15weights_ref.reweight_pkl = gridpack_dir + "18052018_ref/ttgamma/order2/ttgamma0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttgammaLarge_LO_order2_15weights_ref.nEvents      = 10000000
fwlite_ttgammaLarge_LO_order2_15weights_ref.xsec         = 7.092 * 3.697 / 2.302 #pb ttgamma gridpack * ttgamma NLO Daniel / ttgamma LO run.py UFO
#check

#fwlite_ttgammaLarge_LO_order2_15weights_ref_delphes      = Sample.fromDirectory("fwlite_ttgammaLarge_LO_order2_15weights_ref_delphes", texName = "ttgamma",  directory = [os.path.join( gen_dir, "fwlite_ttgammaLarge_LO_order2_15weights_ref", "delphes" )], treeName = "Delphes")
#fwlite_ttgammaLarge_LO_order2_15weights_ref.addFriend( fwlite_ttgammaLarge_LO_order2_15weights_ref_delphes, treeName = "Delphes", sortFiles = True)

##############################
##############################
##############################

# signal samples with reference point

#gen_dir = "/afs/hephy.at/data/llechner01/TTXPheno/skims/gen/v17/"

# reference point 15/2
fwlite_ttZ_ll_LO_order2_15weights_ref               = Sample.fromDirectory("fwlite_ttZ_ll_LO_order2_15weights_ref",  texName = "ttZ",      directory = [os.path.join( gen_dir, "fwlite_ttZ_ll_LO_order2_15weights_ref" )]) 
fwlite_ttZ_ll_LO_order2_15weights_ref.reweight_pkl  = gridpack_dir + "18052018_ref/ttZ/order2/ttZ0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttZ_ll_LO_order2_15weights_ref.nEvents       = 990000
fwlite_ttZ_ll_LO_order2_15weights_ref.xsec          = 0.5205 * 0.0915 / 0.0565 #pb ttZ, Z->ll, ttZ gridpack * ttZ NLO Daniel / ttZ LO run.py UFO
#check

fwlite_ttW_LO_order2_15weights_ref                  = Sample.fromDirectory("fwlite_ttW_LO_order2_15weights_ref",     texName = "ttW",      directory = [os.path.join( gen_dir, "fwlite_ttW_LO_order2_15weights_ref" )])
fwlite_ttW_LO_order2_15weights_ref.reweight_pkl     = gridpack_dir + "18052018_ref/ttW/order2/ttW0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttW_LO_order2_15weights_ref.nEvents          = 945000
fwlite_ttW_LO_order2_15weights_ref.xsec             = 0.3599 * 0.2043 / 0.1336 #pb ttW, ttW gridpack * ttW NLO Daniel / ttW LO run.py UFO *BR
#check

fwlite_ttgamma_LO_order2_15weights_ref              = Sample.fromDirectory("fwlite_ttgamma_LO_order2_15weights_ref", texName = "ttgamma",  directory = [os.path.join( gen_dir, "fwlite_ttgamma_LO_order2_15weights_ref" )])
fwlite_ttgamma_LO_order2_15weights_ref.reweight_pkl = gridpack_dir + "18052018_ref/ttgamma/order2/ttgamma0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttgamma_LO_order2_15weights_ref.nEvents      = 970000
fwlite_ttgamma_LO_order2_15weights_ref.xsec         = 7.092 * 3.697 / 2.302 #pb ttgamma gridpack * ttgamma NLO Daniel / ttgamma LO run.py UFO
#check

#fwlite_ttZ_ll_LO_order2_15weights_ref_delphes       = Sample.fromDirectory("fwlite_ttZ_ll_LO_order2_15weights_ref_delphes",  texName = "ttZ",      directory = [os.path.join( gen_dir, "fwlite_ttZ_ll_LO_order2_15weights_ref", "delphes" )], treeName = "Delphes") 
#fwlite_ttW_LO_order2_15weights_ref_delphes          = Sample.fromDirectory("fwlite_ttW_LO_order2_15weights_ref_delphes",     texName = "ttW",      directory = [os.path.join( gen_dir, "fwlite_ttW_LO_order2_15weights_ref", "delphes" )], treeName = "Delphes")
#fwlite_ttgamma_LO_order2_15weights_ref_delphes      = Sample.fromDirectory("fwlite_ttgamma_LO_order2_15weights_ref_delphes", texName = "ttgamma",  directory = [os.path.join( gen_dir, "fwlite_ttgamma_LO_order2_15weights_ref", "delphes" )], treeName = "Delphes")

#fwlite_ttZ_ll_LO_order2_15weights_ref.addFriend(  fwlite_ttZ_ll_LO_order2_15weights_ref_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_ttW_LO_order2_15weights_ref.addFriend(     fwlite_ttW_LO_order2_15weights_ref_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_ttgamma_LO_order2_15weights_ref.addFriend( fwlite_ttgamma_LO_order2_15weights_ref_delphes, treeName = "Delphes", sortFiles = True)

##############################
##############################
##############################

# backgrounds without reference point
# leptonic decays W > lnu, Z > ll, t > Wb
fwlite_tt_dilep_LO_order2_15weights               = Sample.fromDirectory("fwlite_tt_dilep_LO_order2_15weights",  texName = "ttdilep",      directory = [os.path.join( gen_dir, "fwlite_tt_dilep_LO_order2_15weights" )]) 
fwlite_tt_dilep_LO_order2_15weights.reweight_pkl  = gridpack_dir + "06072018/tt/order2/tt_lep_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_tt_dilep_LO_order2_15weights.nEvents       = 1000000 
fwlite_tt_dilep_LO_order2_15weights.xsec          = 45.71 * 831.76 / 485.8 #tt dilep gridpack * tt NLO Daniel / tt LO run.py UFO
#check

fwlite_tt_nonhad_LO_order2_15weights               = Sample.fromDirectory("fwlite_tt_nonhad_LO_order2_15weights",  texName = "ttnonhad",      directory = [os.path.join( gen_dir, "fwlite_tt_nonhad_LO_order2_15weights" )]) 
fwlite_tt_nonhad_LO_order2_15weights.reweight_pkl  = gridpack_dir + "06072018/tt_semilep/order2/tt_semilep_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_tt_nonhad_LO_order2_15weights.nEvents       = 1000000 
fwlite_tt_nonhad_LO_order2_15weights.xsec          = 300.7 * 831.76 / 485.8 #tt nonhad gridpack * tt NLO Daniel / tt LO run.py UFO
#check

fwlite_tt_full_LO_order2_15weights               = Sample.fromDirectory("fwlite_tt_full_LO_order2_15weights",  texName = "tt",      directory = [os.path.join( gen_dir, "fwlite_tt_full_LO_order2_15weights" )]) 
fwlite_tt_full_LO_order2_15weights.reweight_pkl  = gridpack_dir + "06072018/tt_semilep/order2/tt_semilep_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_tt_full_LO_order2_15weights.nEvents       = 10000000 
fwlite_tt_full_LO_order2_15weights.xsec          = 494.9 * 831.76 / 485.8 #tt full gridpack * tt NLO Daniel / tt LO run.py UFO
#check

#fwlite_tt_dilep_LO_order2_15weights_delphes     = Sample.fromDirectory("fwlite_tt_dilep_LO_order2_15weights_ref_delphes",  texName = "ttdilep",      directory = [os.path.join( gen_dir, "fwlite_tt_dilep_LO_order2_15weights", "delphes" )], treeName = "Delphes") 
#fwlite_tt_nonhad_LO_order2_15weights_delphes    = Sample.fromDirectory("fwlite_tt_nonhad_LO_order2_15weights_ref_delphes", texName = "ttnonhad",     directory = [os.path.join( gen_dir, "fwlite_tt_nonhad_LO_order2_15weights", "delphes" )], treeName = "Delphes") 
#fwlite_tt_full_LO_order2_15weights_delphes      = Sample.fromDirectory("fwlite_tt_full_LO_order2_15weights_ref_delphes",   texName = "tt",           directory = [os.path.join( gen_dir, "fwlite_tt_full_LO_order2_15weights", "delphes" )], treeName = "Delphes") 

#fwlite_tt_full_LO_order2_15weights.addFriend(    fwlite_tt_full_LO_order2_15weights_ref_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_tt_nonhad_LO_order2_15weights.addFriend(  fwlite_tt_nonhad_LO_order2_15weights_ref_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_tt_dilep_LO_order2_15weights.addFriend(   fwlite_tt_dilep_LO_order2_15weights_ref_delphes, treeName = "Delphes", sortFiles = True)

##############################
##############################
##############################

fwlite_tZq_LO_order2_15weights               = Sample.fromDirectory("fwlite_tZq_LO_order2_15weights",  texName = "tZq",      directory = [os.path.join( gen_dir, "fwlite_tZq_LO_order2_15weights" )]) 
fwlite_tZq_LO_order2_15weights.reweight_pkl  = gridpack_dir + "06072018/tZq/order2/tZq_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_tZq_LO_order2_15weights.nEvents       = 1000000 
fwlite_tZq_LO_order2_15weights.xsec          = 0.06774 * 0.0758 / 0.06148 #pb tZq gridpack * tZq NLO Daniel / tZq LO run.py DIM6
#check

fwlite_tW_LO_order2_15weights               = Sample.fromDirectory("fwlite_tW_LO_order2_15weights",  texName = "tW",      directory = [os.path.join( gen_dir, "fwlite_tW_LO_order2_15weights" )]) 
fwlite_tW_LO_order2_15weights.reweight_pkl  = gridpack_dir + "06072018/tW/order2/tW_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_tW_LO_order2_15weights.nEvents       = 1000000 
fwlite_tW_LO_order2_15weights.xsec          = 15.82 * 19.55 / 18.09 #pb tW gridpack * tW NLO Daniel / tW LO run.py UFO
#check

fwlite_tWZ_LO_order2_15weights               = Sample.fromDirectory("fwlite_tWZ_LO_order2_15weights",  texName = "tWZ",      directory = [os.path.join( gen_dir, "fwlite_tWZ_LO_order2_15weights" )]) 
fwlite_tWZ_LO_order2_15weights.reweight_pkl  = gridpack_dir + "06072018/tWZ/order2/tWZ_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_tWZ_LO_order2_15weights.nEvents       = 1000000 
fwlite_tWZ_LO_order2_15weights.xsec          = 0.01035 * 0.01123 / 0.01225 #tWZ gridpack * tWll NLO Daniel / tWll LO run.py UFO
#check

fwlite_Zgamma_LO_order2_15weights               = Sample.fromDirectory("fwlite_Zgamma_LO_order2_15weights",  texName = "Zgamma",      directory = [os.path.join( gen_dir, "fwlite_Zgamma_LO_order2_15weights" )]) 
fwlite_Zgamma_LO_order2_15weights.reweight_pkl  = gridpack_dir + "06072018/Zgamma/order2/Zgamma_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_Zgamma_LO_order2_15weights.nEvents       = 1000000 
fwlite_Zgamma_LO_order2_15weights.xsec          = 131.3 #pb ZGTo2LG_ext Daniel
#check?

fwlite_ttgamma_bg_LO_order2_15weights               = Sample.fromDirectory("fwlite_ttgamma_bg_LO_order2_15weights",  texName = "ttgamma",      directory = [os.path.join( gen_dir, "fwlite_ttgamma_bg_LO_order2_15weights" )]) 
fwlite_ttgamma_bg_LO_order2_15weights.reweight_pkl  = gridpack_dir + "06072018/ttgamma/order2/ttgamma_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttgamma_bg_LO_order2_15weights.nEvents       = 1000000 
fwlite_ttgamma_bg_LO_order2_15weights.xsec          = 2.179 * 3.697 / 2.302 #pb ttgamma gridpack * ttgamma NLO Daniel / ttgamma LO run.py UFO
#check

fwlite_WZ_lep_LO_order2_15weights                   = Sample.fromDirectory("fwlite_WZ_lep_LO_order2_15weights",     texName = "WZ",      directory = [os.path.join( gen_dir, "fwlite_WZ_lep_LO_order2_15weights" )])
fwlite_WZ_lep_LO_order2_15weights.reweight_pkl      = gridpack_dir + "06072018/WZ/order2/WZ_lep_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_WZ_lep_LO_order2_15weights.nEvents           = 1000000 
fwlite_WZ_lep_LO_order2_15weights.xsec              = 4.666 #WZTo3LNu_amcatnlo  #47.13*(3*0.108)*(3*0.0336) #pb WZ NLO Daniel * BR(Wlep) * BR(Zlep)
#check?

#fwlite_tZq_LO_order2_15weights_delphes           = Sample.fromDirectory("fwlite_tZq_LO_order2_15weights_ref_delphes",         texName = "tZq",        directory = [os.path.join( gen_dir, "fwlite_tZq_LO_order2_15weights", "delphes" )], treeName = "Delphes") 
#fwlite_tW_LO_order2_15weights_delphes            = Sample.fromDirectory("fwlite_tW_LO_order2_15weights_ref_delphes",          texName = "tW",         directory = [os.path.join( gen_dir, "fwlite_tW_LO_order2_15weights", "delphes" )], treeName = "Delphes") 
#fwlite_tWZ_LO_order2_15weights_delphes           = Sample.fromDirectory("fwlite_tWZ_LO_order2_15weights_ref_delphes",         texName = "tWZ",        directory = [os.path.join( gen_dir, "fwlite_tWZ_LO_order2_15weights", "delphes" )], treeName = "Delphes") 
#fwlite_Zgamma_LO_order2_15weights_delphes        = Sample.fromDirectory("fwlite_Zgamma_LO_order2_15weights_ref_delphes",      texName = "Zgamma",     directory = [os.path.join( gen_dir, "fwlite_Zgamma_LO_order2_15weights", "delphes" )], treeName = "Delphes") 
#fwlite_ttgamma_bg_LO_order2_15weights_delphes    = Sample.fromDirectory("fwlite_ttgamma_bg_LO_order2_15weights_ref_delphes",  texName = "ttgamma_bg", directory = [os.path.join( gen_dir, "fwlite_ttgamma_bg_LO_order2_15weights", "delphes" )], treeName = "Delphes") 
#fwlite_WZ_lep_LO_order2_15weights_delphes        = Sample.fromDirectory("fwlite_WZ_lep_LO_order2_15weights_delphes",      texName = "WZ",      directory = [os.path.join( gen_dir, "fwlite_WZ_lep_LO_order2_15weights", "delphes" )], treeName = "Delphes")

#fwlite_tZq_LO_order2_15weights.addFriend(        fwlite_tZq_LO_order2_15weights_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_tW_LO_order2_15weights.addFriend(         fwlite_tW_LO_order2_15weights_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_tWZ_LO_order2_15weights.addFriend(        fwlite_tWz_LO_order2_15weights_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_Zgamma_LO_order2_15weights.addFriend(     fwlite_Zgamma_LO_order2_15weights_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_ttgamma_bg_LO_order2_15weights.addFriend( fwlite_ttgamma_bg_LO_order2_15weights_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_WZ_lep_LO_order2_15weights.addFriend(     fwlite_WZ_lep_LO_order2_15weights_delphes, treeName = "Delphes", sortFiles = True)

##############################
##############################
##############################

##############################
##############################
##############################

"""
# Signal samples no reference point (not indended to use)

#gen_dir = "/afs/hephy.at/data/rschoefbeck01/TTXPheno/skims/gen/v18/"
gen_dir = "/afs/hephy.at/data/llechner01/TTXPheno/skims/gen/v17/"

# no reference point 15/2
fwlite_ttZ_ll_LO_order2_15weights               = Sample.fromDirectory("fwlite_ttZ_ll_LO_order2_15weights",  texName = "ttZ",      directory = [os.path.join( gen_dir, "fwlite_ttZ_ll_LO_order2_15weights" )]) 
fwlite_ttZ_ll_LO_order2_15weights.reweight_pkl  = gridpack_dir + "18052018/ttZ/order2/ttZ0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttZ_ll_LO_order2_15weights.nEvents       = 975000
fwlite_ttZ_ll_LO_order2_15weights.xsec          = 0.05363 * 0.0915 / 0.0565 #pb ttZ, Z->ll, ttZ gridpack * ttZ NLO Daniel / ttZ LO run.py UFO
#check

#fwlite_ttW_LO_order2_15weights                  = Sample.fromDirectory("fwlite_ttW_LO_order2_15weights",     texName = "ttW",      directory = [os.path.join( gen_dir, "fwlite_ttW_LO_order2_15weights" )])
#fwlite_ttW_LO_order2_15weights.reweight_pkl     = gridpack_dir + "18052018/ttW/order2/ttW0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
#fwlite_ttW_LO_order2_15weights.nEvents          = 1000000
#fwlite_ttW_LO_order2_15weights.xsec             = 0.1134 * 0.2043 / 0.1336 #pb ttW, ttW gridpack * ttW NLO Daniel / ttW LO run.py UFO *BR
#check

fwlite_ttgamma_LO_order2_15weights              = Sample.fromDirectory("fwlite_ttgamma_LO_order2_15weights", texName = "ttgamma",  directory = [os.path.join( gen_dir, "fwlite_ttgamma_LO_order2_15weights" )])
fwlite_ttgamma_LO_order2_15weights.reweight_pkl = gridpack_dir + "18052018/ttgamma/order2/ttgamma0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttgamma_LO_order2_15weights.nEvents      = 990000
fwlite_ttgamma_LO_order2_15weights.xsec         = 2.176 * 3.697 / 2.302 #pb ttgamma gridpack * ttZ NLO Daniel / ttgamma LO run.py UFO
#check

fwlite_ttZ_ll_LO_order2_15weights_delphes       = Sample.fromDirectory("fwlite_ttZ_ll_LO_order2_15weights_delphes",  texName = "ttZ",      directory = [os.path.join( gen_dir, "fwlite_ttZ_ll_LO_order2_15weights", "delphes" )], treeName = "Delphes") 
#fwlite_ttW_LO_order2_15weights_delphes          = Sample.fromDirectory("fwlite_ttW_LO_order2_15weights_delphes",     texName = "ttW",      directory = [os.path.join( gen_dir, "fwlite_ttW_LO_order2_15weights", "delphes" )], treeName = "Delphes")
fwlite_ttgamma_LO_order2_15weights_delphes      = Sample.fromDirectory("fwlite_ttgamma_LO_order2_15weights_delphes", texName = "ttgamma",  directory = [os.path.join( gen_dir, "fwlite_ttgamma_LO_order2_15weights", "delphes" )], treeName = "Delphes")

#fwlite_ttZ_ll_LO_order2_15weights.addFriend(  fwlite_ttZ_ll_LO_order2_15weights_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_ttW_LO_order2_15weights.addFriend(     fwlite_ttW_LO_order2_15weights_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_ttgamma_LO_order2_15weights.addFriend( fwlite_ttgamma_LO_order2_15weights_delphes, treeName = "Delphes", sortFiles = True)

##############################
##############################
##############################

# signal sample without reference point with 8 weights, 3rd order (not indended to use)

gen_dir = "/afs/hephy.at/data/llechner01/TTXPheno/skims/gen/v3/"

# no reference point 8/3
fwlite_ttZ_ll_LO_order3_8weights               = Sample.fromDirectory("fwlite_ttZ_ll_LO_order3_8weights",  texName = "ttZ",      directory = [os.path.join( gen_dir, "fwlite_ttZ_ll_LO_order3_8weights" )]) 
fwlite_ttZ_ll_LO_order3_8weights.reweight_pkl  = gridpack_dir + "07052018/ttZ/order3/ttZ0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttZ_ll_LO_order3_8weights.nEvents       = 995000
fwlite_ttZ_ll_LO_order3_8weights.xsec          = 0.05363 * 0.0915 / 0.0565 #pb ttZ, Z->ll, ttZ gridpack * ttZ NLO Daniel / ttZ LO run.py UFO
#check

fwlite_ttW_LO_order3_8weights                  = Sample.fromDirectory("fwlite_ttW_LO_order3_8weights",     texName = "ttW",      directory = [os.path.join( gen_dir, "fwlite_ttW_LO_order3_8weights" )])
fwlite_ttW_LO_order3_8weights.reweight_pkl     = gridpack_dir + "07052018/ttW/order3/ttW0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttW_LO_order3_8weights.nEvents          = 1000000
fwlite_ttW_LO_order3_8weights.xsec             = 0.1134 * 0.2043 / 0.1336 #pb ttW, ttW gridpack * ttW NLO Daniel / ttW LO run.py UFO *BR
#check

fwlite_ttgamma_LO_order3_8weights              = Sample.fromDirectory("fwlite_ttgamma_LO_order3_8weights", texName = "ttgamma",  directory = [os.path.join( gen_dir, "fwlite_ttgamma_LO_order3_8weights" )])
fwlite_ttgamma_LO_order3_8weights.reweight_pkl = gridpack_dir + "07052018/ttgamma/order3/ttgamma0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttgamma_LO_order3_8weights.nEvents      = 1000000
fwlite_ttgamma_LO_order3_8weights.xsec         = 2.176 * 3.697 / 2.302 #pb ttgamma gridpack * ttZ NLO Daniel / ttgamma LO run.py UFO
#check

fwlite_ttZ_ll_LO_order3_8weights_delphes       = Sample.fromDirectory("fwlite_ttZ_ll_LO_order3_8weights_delphes",  texName = "ttZ",      directory = [os.path.join( gen_dir, "fwlite_ttZ_ll_LO_order3_8weights", "delphes" )], treeName = "Delphes") 
fwlite_ttW_LO_order3_8weights_delphes          = Sample.fromDirectory("fwlite_ttW_LO_order3_8weights_delphes",     texName = "ttW",      directory = [os.path.join( gen_dir, "fwlite_ttW_LO_order3_8weights", "delphes" )], treeName = "Delphes")
fwlite_ttgamma_LO_order3_8weights_delphes      = Sample.fromDirectory("fwlite_ttgamma_LO_order3_8weights_delphes", texName = "ttgamma",  directory = [os.path.join( gen_dir, "fwlite_ttgamma_LO_order3_8weights", "delphes" )], treeName = "Delphes")

#fwlite_ttZ_ll_LO_order3_8weights.addFriend(  fwlite_ttZ_ll_LO_order3_8weights_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_ttW_LO_order3_8weights.addFriend(     fwlite_ttW_LO_order3_8weights_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_ttgamma_LO_order3_8weights.addFriend( fwlite_ttgamma_LO_order3_8weights_delphes, treeName = "Delphes", sortFiles = True)

##############################
##############################
##############################

# full ttbar sample with reference point (not indended to use)
# hadronic + leptonic decays
fwlite_tt_LO_order2_15weights_ref               = Sample.fromDirectory("fwlite_tt_LO_order2_15weights_ref",  texName = "tt",      directory = [os.path.join( gen_dir, "fwlite_tt_LO_order2_15weights_ref" )]) 
fwlite_tt_LO_order2_15weights_ref.reweight_pkl  = gridpack_dir + "04062018/tt/order2/tt_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_tt_LO_order2_15weights_ref.nEvents       = 1000000 
fwlite_tt_LO_order2_15weights_ref.xsec          = 1591 * 831.76 / 485.8 #pb tt, tt dilep xsec BSM LO * tt NLO Daniel / tt LO run.py UFO
#check

fwlite_tt_LO_order2_15weights_ref_delphes       = Sample.fromDirectory("fwlite_tt_LO_order2_15weights_ref_delphes",  texName = "tt",      directory = [os.path.join( gen_dir, "fwlite_tt_LO_order2_15weights_ref", "delphes" )], treeName = "Delphes") 

# dilep ttbar sample with reference point (not indended to use)
# leptonic decays W > lnu, Z > ll, t > Wb
fwlite_tt_dilep_LO_order2_15weights_ref               = Sample.fromDirectory("fwlite_tt_dilep_LO_order2_15weights_ref",  texName = "tt",      directory = [os.path.join( gen_dir, "fwlite_tt_dilep_LO_order2_15weights_ref" )]) 
fwlite_tt_dilep_LO_order2_15weights_ref.reweight_pkl  = gridpack_dir + "05062018_ref/tt/order2/tt_lep_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_tt_dilep_LO_order2_15weights_ref.nEvents       = 1000000 
fwlite_tt_dilep_LO_order2_15weights_ref.xsec          = 174.6 * 831.76 / 485.8 #pb tt, tt dilep xsec BSM LO * tt NLO Daniel / tt LO run.py UFO
#check

# full WZ sample (not indended to use)
fwlite_WZ_LO_order2_15weights                   = Sample.fromDirectory("fwlite_WZ_LO_order2_15weights",     texName = "WZ",      directory = [os.path.join( gen_dir, "fwlite_WZ_LO_order2_15weights" )])
fwlite_WZ_LO_order2_15weights.reweight_pkl      = gridpack_dir + "04062018/WZ/order2/WZ_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_WZ_LO_order2_15weights.nEvents           = 1000000 
fwlite_WZ_LO_order2_15weights.xsec              = 47.13 #pb WZ, xsec NLO Daniel
#check?

fwlite_tt_dilep_LO_order2_15weights_ref_delphes   = Sample.fromDirectory("fwlite_tt_dilep_LO_order2_15weights_ref_delphes",  texName = "tt",      directory = [os.path.join( gen_dir, "fwlite_tt_dilep_LO_order2_15weights_ref", "delphes" )], treeName = "Delphes") 
fwlite_WZ_LO_order2_15weights_delphes           = Sample.fromDirectory("fwlite_WZ_LO_order2_15weights_delphes",      texName = "WZ",      directory = [os.path.join( gen_dir, "fwlite_WZ_LO_order2_15weights", "delphes" )], treeName = "Delphes")

#fwlite_tt_dilep_LO_order2_15weights_ref.addFriend(  fwlite_tt_dilep_LO_order2_15weights_ref_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_tt_LO_order2_15weights_ref.addFriend(  fwlite_tt_LO_order2_15weights_ref_delphes, treeName = "Delphes", sortFiles = True)
#fwlite_WZ_LO_order2_15weights.addFriend(      fwlite_WZ_LO_order2_15weights_delphes, treeName = "Delphes", sortFiles = True)

##############################
##############################
##############################

# old test samples
gen_dir = "/afs/hephy.at/data/rschoefbeck02/TopEFT/skims/gen/v2/"

dim6top_ttZ_ll_LO_highStat_scan                           = Sample.fromDirectory("fwlite_ttZ_ll_LO_highStat_scan", texName = "ttZ (scan)", directory = [os.path.join( gen_dir, "fwlite_ttZ_ll_LO_highStat_scan")])
dim6top_ttZ_ll_LO_highStat_scan.reweight_pkl              = "/afs/hephy.at/data/rschoefbeck02/TopEFT/results/gridpacks/ttZ0j_rwgt_patch_625_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"

dim6top_ttZ_ll_LO_currentplane_highStat_scan              = Sample.fromDirectory("dim6top_ttZ_ll_LO_currentplane_highStat_scan", texName = "ttZ (current scan)", directory = [os.path.join( gen_dir, "fwlite_ttZ_ll_LO_currentplane_highStat_scan")])
dim6top_ttZ_ll_LO_currentplane_highStat_scan.reweight_pkl = "/afs/hephy.at/data/rschoefbeck02/TopEFT/results/gridpacks/ttZ0j_rwgt_patch_currentplane_highStat_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"

test = Sample.fromFiles("test", files = ["/afs/hephy.at/data/rschoefbeck02/TTXPheno/skims/gen/v2/test/test.root"], texName = "test")
test.reweight_pkl = '/afs/cern.ch/user/l/llechner/public/gridpacks_data/order_3/ttZ/gridpacks/reweight_card.pkl'

"""
