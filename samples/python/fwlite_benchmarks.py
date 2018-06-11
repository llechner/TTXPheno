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
test.nEvents = 10
test.xsec = 200


# no reference point samples 8/3
fwlite_ttZ_ll_LO_order3_8weights               = FWLiteSample.fromDAS("fwlite_ttZ_ll_LO_order3_8weights", "/ttZ0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball/llechner-ttZ0j_order3_8weights-7a5fde3f5bf89006ee3acec926ca87d8/USER", "phys03", dbFile = dbFile)
fwlite_ttZ_ll_LO_order3_8weights.reweight_pkl  = "/afs/hephy.at/data/llechner01/TTXPheno/gridpacks/07052018/ttZ/order3/ttZ0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttZ_ll_LO_order3_8weights.nEvents       = 995000
fwlite_ttZ_ll_LO_order3_8weights.xsec          = 0.0915 #pb ttZ, Z->ll

fwlite_ttW_LO_order3_8weights                  = FWLiteSample.fromDAS("fwlite_ttW_LO_order3_8weights", "/ttW0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball/llechner-ttW0j_order3_8weights-593ea75549b4c51667dffc93040bbda1/USER", "phys03", dbFile = dbFile)
fwlite_ttW_LO_order3_8weights.reweight_pkl     = "/afs/hephy.at/data/llechner01/TTXPheno/gridpacks/07052018/ttW/order3/ttW0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttW_LO_order3_8weights.nEvents          = 1000000
fwlite_ttW_LO_order3_8weights.xsec             = 0.2043 #pb ttW, W->lnu

fwlite_ttgamma_LO_order3_8weights              = FWLiteSample.fromDAS("fwlite_ttgamma_LO_order3_8weights", "/ttgamma0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball/llechner-ttgamma0j_order3_8weights-10fcfa1a1c01204983ea66975abf2caf/USER", "phys03", dbFile = dbFile)
fwlite_ttgamma_LO_order3_8weights.reweight_pkl = "/afs/hephy.at/data/llechner01/TTXPheno/gridpacks/07052018/ttgamma/order3/ttgamma0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttgamma_LO_order3_8weights.nEvents      = 1000000
fwlite_ttgamma_LO_order3_8weights.xsec         = 3.697 #pb ttgamma


# no reference point samples 15/2
fwlite_ttZ_ll_LO_order2_15weights               = FWLiteSample.fromDAS("fwlite_ttZ_ll_LO_order2_15weights", "/ttZ0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball/llechner-ttZ0j_order2_15weights_18052018-7a5fde3f5bf89006ee3acec926ca87d8/USER", "phys03", dbFile = dbFile)
fwlite_ttZ_ll_LO_order2_15weights.reweight_pkl  = "/afs/hephy.at/data/llechner01/TTXPheno/gridpacks/18052018/ttZ/order2/ttZ0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttZ_ll_LO_order2_15weights.nEvents	= 975000
fwlite_ttZ_ll_LO_order2_15weights.xsec          = 0.0915 #pb ttZ, Z->ll

fwlite_ttW_LO_order2_15weights                  = FWLiteSample.fromDAS("fwlite_ttW_LO_order2_15weights", "/ttW0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball/llechner-ttW0j_order2_15weights_18052018-593ea75549b4c51667dffc93040bbda1/USER", "phys03", dbFile = dbFile)
fwlite_ttW_LO_order2_15weights.reweight_pkl     = "/afs/hephy.at/data/llechner01/TTXPheno/gridpacks/18052018/ttW/order2/ttW0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttW_LO_order2_15weights.nEvents          = 1000000
fwlite_ttW_LO_order2_15weights.xsec             = 0.2043 #pb ttW, W->lnu

fwlite_ttgamma_LO_order2_15weights              = FWLiteSample.fromDAS("fwlite_ttgamma_LO_order2_15weights", "/ttgamma0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball/llechner-ttgamma0j_order2_15weights_18052018-10fcfa1a1c01204983ea66975abf2caf/USER", "phys03", dbFile = dbFile)
fwlite_ttgamma_LO_order2_15weights.reweight_pkl = "/afs/hephy.at/data/llechner01/TTXPheno/gridpacks/18052018/ttgamma/order2/ttgamma0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttgamma_LO_order2_15weights.nEvents	= 990000
fwlite_ttgamma_LO_order2_15weights.xsec         = 3.697 #pb ttgamma


# reference point samples 15/2
fwlite_ttZ_ll_LO_order2_15weights_ref               = FWLiteSample.fromDAS("fwlite_ttZ_ll_LO_order2_15weights_ref", "/ttZ0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball/llechner-ttZ0j_order2_15weights_18052018_ref-7a5fde3f5bf89006ee3acec926ca87d8/USER", "phys03", dbFile = dbFile)
fwlite_ttZ_ll_LO_order2_15weights_ref.reweight_pkl  = "/afs/hephy.at/data/llechner01/TTXPheno/gridpacks/18052018_ref/ttZ/order2/ttZ0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttZ_ll_LO_order2_15weights_ref.nEvents       = 990000
fwlite_ttZ_ll_LO_order2_15weights_ref.xsec          = 0.8815 #pb ttZ, Z->ll, xsec_SM_NNLO * xsec_BSM_LO / xsec_SM_LO

fwlite_ttW_LO_order2_15weights_ref                  = FWLiteSample.fromDAS("fwlite_ttW_LO_order2_15weights_ref", "/ttW0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball/llechner-ttW0j_order2_15weights_18052018_ref-593ea75549b4c51667dffc93040bbda1/USER", "phys03", dbFile = dbFile)
fwlite_ttW_LO_order2_15weights_ref.reweight_pkl     = "/afs/hephy.at/data/llechner01/TTXPheno/gridpacks/18052018_ref/ttW/order2/ttW0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttW_LO_order2_15weights_ref.nEvents          = 945000
fwlite_ttW_LO_order2_15weights_ref.xsec             = 0.6327 #pb ttW, W->lnu, xsec_SM_NNLO * xsec_BSM_LO / xsec_SM_LO

fwlite_ttgamma_LO_order2_15weights_ref              = FWLiteSample.fromDAS("fwlite_ttgamma_LO_order2_15weights_ref", "/ttgamma0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball/llechner-ttgamma0j_order2_15weights_18052018_ref-10fcfa1a1c01204983ea66975abf2caf/USER", "phys03", dbFile = dbFile)
fwlite_ttgamma_LO_order2_15weights_ref.reweight_pkl = "/afs/hephy.at/data/llechner01/TTXPheno/gridpacks/18052018_ref/ttgamma/order2/ttgamma0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_ttgamma_LO_order2_15weights_ref.nEvents      = 970000
fwlite_ttgamma_LO_order2_15weights_ref.xsec         = 11.8807 #pb ttgamma, xsec_SM_NNLO * xsec_BSM_LO / xsec_SM_LO

# background samples 0/2 (WZ) or 15/2 (ttbar with reference point)
#leptonic decays W > lnu, Z > ll, t > Wb
fwlite_tt_lep_LO_order2_15weights_ref              = FWLiteSample.fromDAS("fwlite_tt_lep_LO_order2_15weights_ref", "/tt_lep_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball/llechner-bg_lep_dim6top_05June18-399ed716eb7225402bb4416ff36fe4d6/USER", "phys03", dbFile = dbFile)
fwlite_tt_lep_LO_order2_15weights_ref.reweight_pkl = "/afs/hephy.at/data/llechner01/TTXPheno/gridpacks/05062018/tt/order2/tt_lep_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_tt_lep_LO_order2_15weights_ref.nEvents      = 1000000 #? not checked!
fwlite_tt_lep_LO_order2_15weights_ref.xsec         = 87.3148 * 735.1 / 46.57 #pb tt, W->lnu, xsec_SM_NNLO * xsec_BSM_LO / xsec_SM_LO

fwlite_WZ_lep_LO_order2_15weights                  = FWLiteSample.fromDAS("fwlite_WZ_lep_LO_order2_15weights", "/WZ_lep_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball/llechner-bg_lep_dim6top_05June18-6dade6042d6868e7bb87de85663e8a54/USER", "phys03", dbFile = dbFile)
fwlite_WZ_lep_LO_order2_15weights.reweight_pkl     = "/afs/hephy.at/data/llechner01/TTXPheno/gridpacks/05062018/WZ/order2/WZ_lep_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_WZ_lep_LO_order2_15weights.nEvents          = 1000000 #? not checked!
fwlite_WZ_lep_LO_order2_15weights.xsec             = 4.666 #pb WZ, W->lnu, Z->ll, xsec_SM_NNLO

#hadronic and leptonic decays
fwlite_tt_LO_order2_15weights_ref              = FWLiteSample.fromDAS("fwlite_tt_LO_order2_15weights_ref", "/tt_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball/llechner-bg_dim6top_04June18-db4155f01c4c21dc10125760597b536e/USER", "phys03", dbFile = dbFile)
fwlite_tt_LO_order2_15weights_ref.reweight_pkl = "/afs/hephy.at/data/llechner01/TTXPheno/gridpacks/04062018/tt/order2/tt_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_tt_LO_order2_15weights_ref.nEvents      = 1000000 #? not checked!
fwlite_tt_LO_order2_15weights_ref.xsec         = 2678.803 #pb tt, xsec_SM_NNLO * xsec_BSM_LO / xsec_SM_LO

fwlite_WZ_LO_order2_15weights                  = FWLiteSample.fromDAS("fwlite_WZ_LO_order2_15weights", "/WZ_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball/llechner-bg_dim6top_04June18-a4eb1cf113188459e4c517c4aa52c31c/USER", "phys03", dbFile = dbFile)
fwlite_WZ_LO_order2_15weights.reweight_pkl     = "/afs/hephy.at/data/llechner01/TTXPheno/gridpacks/04062018/WZ/order2/WZ_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
fwlite_WZ_LO_order2_15weights.nEvents          = 1000000 #? not checked!
fwlite_WZ_LO_order2_15weights.xsec             = 47.13 #pb WZ, xsec_SM_NNLO


