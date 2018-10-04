#!/bin/sh

#use cardfile delphes_card_ATLAS or delphes_card_CMS

#########################
##4th order 2WC samples##
#########################

#python genPostProcessing.py --overwrite target --targetDir v18 --logLevel DEBUG --sample fwlite_ttZ_ll_LO_order4_4weights_2WC_ref --addReweights --interpolationOrder 4 --delphes --delphesCard delphes_card_CMS #SPLIT200
#python genPostProcessing.py --overwrite target --targetDir v18 --logLevel DEBUG --sample fwlite_ttgamma_LO_order4_4weights_2WC_ref --addReweights --interpolationOrder 4 --delphes --delphesCard delphes_card_CMS #SPLIT200

###############
##ewkDMGZ_UFO##
###############

#python genPostProcessing.py --overwrite target --targetDir v18 --logLevel DEBUG --sample fwlite_ttZ_ll_LO_order4_4weights_ewkDMGZ_ref --addReweights --interpolationOrder 4 --delphes --delphesCard delphes_card_CMS #SPLIT200
#python genPostProcessing.py --overwrite target --targetDir v18 --logLevel DEBUG --sample fwlite_ttgamma_LO_order4_2weights_ewkDMGZ_ref --addReweights --interpolationOrder 4 --delphes --delphesCard delphes_card_CMS #SPLIT200
#python genPostProcessing.py --overwrite target --targetDir v18 --logLevel DEBUG --sample fwlite_ttZ_ll_LO_order4_4weights_ewkDMGZ_ref --addReweights --interpolationOrder 4 --delphes --delphesCard delphes_card_ATLAS #SPLIT200
#python genPostProcessing.py --overwrite target --targetDir v18 --logLevel DEBUG --sample fwlite_ttgamma_LO_order4_2weights_ewkDMGZ_ref --addReweights --interpolationOrder 4 --delphes --delphesCard delphes_card_ATLAS #SPLIT200

#########
##ATLAS##
#########

# Signal samples with reference point
#python genPostProcessing.py --overwrite target --targetDir v18 --logLevel DEBUG --sample fwlite_ttZ_ll_LO_order2_15weights_ref --addReweights --interpolationOrder 2 --delphes --delphesCard delphes_card_ATLAS #SPLIT200

# Large samples for ttbar and ttgamma
#python genPostProcessing.py --overwrite target --targetDir v18 --logLevel DEBUG --sample fwlite_ttgammaLarge_LO_order2_15weights_ref --addReweights --interpolationOrder 2 --delphes --delphesCard delphes_card_ATLAS #SPLIT2000
#python genPostProcessing.py --overwrite target --targetDir v18 --logLevel DEBUG --sample fwlite_tt_full_LO_order2_15weights --interpolationOrder 2 --delphes --delphesCard delphes_card_ATLAS #SPLIT2000

# background samples without reference point
#python genPostProcessing.py --overwrite target --targetDir v18 --logLevel DEBUG --sample fwlite_tt_dilep_LO_order2_15weights --interpolationOrder 2 --delphes --delphesCard delphes_card_ATLAS #SPLIT200
#python genPostProcessing.py --overwrite target --targetDir v18 --logLevel DEBUG --sample fwlite_tt_nonhad_LO_order2_15weights --interpolationOrder 2 --delphes --delphesCard delphes_card_ATLAS #SPLIT200
#python genPostProcessing.py --overwrite target --targetDir v18 --logLevel DEBUG --sample fwlite_tZq_LO_order2_15weights --interpolationOrder 2 --delphes --delphesCard delphes_card_ATLAS #SPLIT200
#python genPostProcessing.py --overwrite target --targetDir v18 --logLevel DEBUG --sample fwlite_tWZ_LO_order2_15weights --interpolationOrder 2 --delphes --delphesCard delphes_card_ATLAS #SPLIT200
#python genPostProcessing.py --overwrite target --targetDir v18 --logLevel DEBUG --sample fwlite_tW_LO_order2_15weights --interpolationOrder 2 --delphes --delphesCard delphes_card_ATLAS #SPLIT200
#python genPostProcessing.py --overwrite target --targetDir v18 --logLevel DEBUG --sample fwlite_Zgamma_LO_order2_15weights --interpolationOrder 2 --delphes --delphesCard delphes_card_ATLAS #SPLIT200
#python genPostProcessing.py --overwrite target --targetDir v18 --logLevel DEBUG --sample fwlite_ttgamma_bg_LO_order2_15weights --interpolationOrder 2 --delphes --delphesCard delphes_card_ATLAS #SPLIT200
#python genPostProcessing.py --overwrite target --targetDir v18 --logLevel DEBUG --sample fwlite_WZ_lep_LO_order2_15weights --interpolationOrder 2 --delphes --delphesCard delphes_card_ATLAS #SPLIT200

# Signal samples with reference point
#python genPostProcessing.py --overwrite target --targetDir v18 --logLevel DEBUG --sample fwlite_ttgamma_LO_order2_15weights_ref --addReweights --interpolationOrder 2 --delphes --delphesCard delphes_card_ATLAS #SPLIT200
#python genPostProcessing.py --overwrite target --targetDir v18 --logLevel DEBUG --sample fwlite_ttW_LO_order2_15weights_ref --addReweights --interpolationOrder 2 --delphes --delphesCard delphes_card_ATLAS #SPLIT200

#######
##CMS##
#######


# Large samples for ttbar and ttgamma

# Signal samples with reference point
python genPostProcessing.py --overwrite target --targetDir v20 --logLevel DEBUG --sample fwlite_ttZ_ll_LO_order2_15weights_ref --addReweights --interpolationOrder 2  --delphes --delphesCard CMS_PhaseII/CMS_PhaseII_200PU_v03 #SPLIT200

## background samples without reference point
python genPostProcessing.py --overwrite target --targetDir v20 --logLevel DEBUG --sample fwlite_WZ_lep_LO_order2_15weights --interpolationOrder 2  --delphes --delphesCard CMS_PhaseII/CMS_PhaseII_200PU_v03 #SPLIT200
python genPostProcessing.py --overwrite target --targetDir v20 --logLevel DEBUG --sample fwlite_tt_dilep_LO_order2_15weights --interpolationOrder 2  --delphes --delphesCard CMS_PhaseII/CMS_PhaseII_200PU_v03 #SPLIT200
python genPostProcessing.py --overwrite target --targetDir v20 --logLevel DEBUG --sample fwlite_tt_nonhad_LO_order2_15weights --interpolationOrder 2  --delphes --delphesCard CMS_PhaseII/CMS_PhaseII_200PU_v03 #SPLIT200
python genPostProcessing.py --overwrite target --targetDir v20 --logLevel DEBUG --sample fwlite_tZq_LO_order2_15weights --interpolationOrder 2  --delphes --delphesCard CMS_PhaseII/CMS_PhaseII_200PU_v03 #SPLIT200
python genPostProcessing.py --overwrite target --targetDir v20 --logLevel DEBUG --sample fwlite_tWZ_LO_order2_15weights --interpolationOrder 2  --delphes --delphesCard CMS_PhaseII/CMS_PhaseII_200PU_v03 #SPLIT200
python genPostProcessing.py --overwrite target --targetDir v20 --logLevel DEBUG --sample fwlite_tW_LO_order2_15weights --interpolationOrder 2  --delphes --delphesCard CMS_PhaseII/CMS_PhaseII_200PU_v03 #SPLIT200

## ttgamma
##python genPostProcessing.py --overwrite target --targetDir v18 --logLevel DEBUG --sample fwlite_Zgamma_LO_order2_15weights --interpolationOrder 2  --delphes --delphesCard CMS_PhaseII/CMS_PhaseII_200PU_v03 #SPLIT200
##python genPostProcessing.py --overwrite target --targetDir v18 --logLevel DEBUG --sample fwlite_ttgamma_bg_LO_order2_15weights --interpolationOrder 2  --delphes --delphesCard CMS_PhaseII/CMS_PhaseII_200PU_v03 #SPLIT200
#python genPostProcessing.py --overwrite target --targetDir v18 --logLevel DEBUG --sample fwlite_ttgammaLarge_LO_order2_15weights_ref --addReweights --interpolationOrder 2  --delphes --delphesCard CMS_PhaseII/CMS_PhaseII_200PU_v03 #SPLIT2000
#python genPostProcessing.py --overwrite target --targetDir v18 --logLevel DEBUG --sample fwlite_tt_full_LO_order2_15weights --interpolationOrder 2  --delphes --delphesCard CMS_PhaseII/CMS_PhaseII_200PU_v03 #SPLIT2000

## Signal samples with reference point
#python genPostProcessing.py --overwrite target --targetDir v18 --removeDelphesFiles --logLevel DEBUG --sample fwlite_ttgamma_LO_order2_15weights_ref --addReweights --interpolationOrder 2  --delphesCard delphes_card_CMS #SPLIT200
#python genPostProcessing.py --overwrite target --targetDir v18 --removeDelphesFiles --logLevel DEBUG --sample fwlite_ttW_LO_order2_15weights_ref --addReweights --interpolationOrder 2  --delphesCard delphes_card_CMS #SPLIT200

######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################

# ttbar sample with reference point (not intended to use)
#python genPostProcessing.py --overwrite target --targetDir v18 --logLevel DEBUG --sample fwlite_tt_LO_order2_15weights_ref --addReweights --interpolationOrder 2 --delphes #SPLIT200

# Signal sample without reference point (not intended to use)
#python genPostProcessing.py --overwrite target --targetDir v18 --logLevel DEBUG --sample fwlite_ttZ_ll_LO_order2_15weights --interpolationOrder 2 --delphes #SPLIT200
#python genPostProcessing.py --overwrite target --targetDir v18 --logLevel DEBUG --sample fwlite_ttgamma_LO_order2_15weights --interpolationOrder 2 --delphes #SPLIT200
#python genPostProcessing.py --overwrite --targetDir v6 --logLevel DEBUG --sample fwlite_ttW_LO_order2_15weights --addReweights --interpolationOrder 2 --delphes #SPLIT200


# Signal sample without reference point, 8 weights, 3rd order (not intended to use)
#python genPostProcessing.py --logLevel DEBUG --overwrite --sample fwlite_ttZ_ll_LO_order3_8weights --addReweights --interpolationOrder 3 --delphes #SPLIT200
#python genPostProcessing.py --logLevel DEBUG --overwrite --sample fwlite_ttW_LO_order3_8weights --addReweights --interpolationOrder 3 --delphes #SPLIT200
#python genPostProcessing.py --logLevel DEBUG --overwrite --sample fwlite_ttgamma_LO_order3_8weights --addReweights --interpolationOrder 3 --delphes #SPLIT200

# Old or test scans
#python genPostProcessing.py --overwrite --sample fwlite_ttZ_ll_LO_highStat_scan --addReweights --delphes #SPLIT200
#python genPostProcessing.py --overwrite --sample fwlite_ttZ_ll_LO_currentplane_highStat_scan --overwrite --addReweights #SPLIT200
#python genPostProcessing.py --sample fwlite_ttZ_ll_LO_scan --overwrite --addReweights #SPLIT200

#python genPostProcessing.py --sample fwlite_ttZ_ll_LO_C2VA_0p2
#python genPostProcessing.py --sample fwlite_ttZ_ll_LO_sm
#
#python genPostProcessing.py --sample fwlite_ttZ_ll_LO_C2VA_0p2
#python genPostProcessing.py --sample fwlite_ttZ_ll_LO_DC1V_0p5_DC1A_0p5
#python genPostProcessing.py --sample fwlite_ttZ_ll_LO_DC1V_m1_DC1A0p5
#
#python genPostProcessing.py --sample fwlite_ttZ_ll_LO_minXSecC1VA_0
#python genPostProcessing.py --sample fwlite_ttZ_ll_LO_minXSecC1VA_1
#python genPostProcessing.py --sample fwlite_ttZ_ll_LO_minXSecC1VA_2
#python genPostProcessing.py --sample fwlite_ttZ_ll_LO_minXSecC1VA_3
#python genPostProcessing.py --sample fwlite_ttZ_ll_LO_minXSecC1VA_4
#python genPostProcessing.py --sample fwlite_ttZ_ll_LO_minXSecC1VA_5
#python genPostProcessing.py --sample fwlite_ttZ_ll_LO_minXSecC1VA_6
#python genPostProcessing.py --sample fwlite_ttZ_ll_LO_minXSecC1VA_7

#python genPostProcessing.py --sample daniel_fwlite_ttZ_ll_LO_sm
#python genPostProcessing.py --sample daniel_fwlite_ttZ_ll_LO_DC2A_0p2_DC2V_0p2
#python genPostProcessing.py --sample daniel_fwlite_ttZ_ll_LO_DC2A_0p4_DC2V_0p4
#python genPostProcessing.py --sample daniel_fwlite_ttZ_ll_LO_DC2A_0p6_DC2V_0p6
#python genPostProcessing.py --sample daniel_fwlite_ttZ_ll_LO_DC2A_0p216_DC2V_0p523
#python genPostProcessing.py --sample daniel_fwlite_ttZ_ll_LO_DC2A_0p523_DC2V_0p216
#python genPostProcessing.py --sample daniel_fwlite_ttZ_ll_LO_DC2A_0p57_DC2V_0p0
#python genPostProcessing.py --sample daniel_fwlite_ttZ_ll_LO_DC2A_0p0_DC2V_0p57

#python genPostProcessing.py --sample daniel_fwlite_ttZ_ll_LO_DC1V_m0p24_DC2V_0p3
#python genPostProcessing.py --sample daniel_fwlite_ttZ_ll_LO_DC1V_m0p24_DC2A_0p28
#python genPostProcessing.py --sample daniel_fwlite_ttZ_ll_LO_DC1A_0p6_DC2V_0p3
#python genPostProcessing.py --sample daniel_fwlite_ttZ_ll_LO_DC1A_0p6_DC2A_0p28

#python genPostProcessing.py --sample daniel_fwlite_ttZ_ll_LO_antiSM_DC2A_0p5
#python genPostProcessing.py --sample daniel_fwlite_ttZ_ll_LO_antiSM_DC2A_0p7

#python genPostProcessing.py --sample daniel_fwlite_ttZ_ll_LO_DC1A_1p19_DC1V_m0p31
#python genPostProcessing.py --sample daniel_fwlite_ttZ_ll_LO_DC1A_0p8_DC1V_0p8
#python genPostProcessing.py --sample daniel_fwlite_ttZ_ll_LO_DC1A_0p6_DC1V_0p59
#python genPostProcessing.py --sample daniel_fwlite_ttZ_ll_LO_DC1A_0p6_DC1V_m1p05
#python genPostProcessing.py --sample daniel_fwlite_ttZ_ll_LO_DC1A_0p6_DC1V_m0p24
#python genPostProcessing.py --sample daniel_fwlite_ttZ_ll_LO_DC1A_0p4_DC1V_0p4
#python genPostProcessing.py --sample daniel_fwlite_ttZ_ll_LO_DC1A_1p0_DC1V_1p0

#python genPostProcessing.py --sample ewkDM_ttZ_ll_GS #SPLIT20
#python genPostProcessing.py --sample ewkDM_ttZ_ll_GS_DC1A_0p500000_DC1V_0p500000 #SPLIT20
#python genPostProcessing.py --sample ewkDM_ttZ_ll_GS_DC1A_0p500000_DC1V_m1p000000 #SPLIT20
#python genPostProcessing.py --sample ewkDM_ttZ_ll_GS_DC1A_0p600000_DC1V_m0p240000_DC2A_0p250000 #SPLIT20
#python genPostProcessing.py --sample ewkDM_ttZ_ll_GS_DC1A_0p600000_DC1V_m0p240000_DC2A_m0p250000 #SPLIT20
#python genPostProcessing.py --sample ewkDM_ttZ_ll_GS_DC1A_0p600000_DC1V_m0p240000_DC2V_0p250000 #SPLIT20
#python genPostProcessing.py --sample ewkDM_ttZ_ll_GS_DC1A_0p600000_DC1V_m0p240000_DC2V_m0p250000 #SPLIT20
#python genPostProcessing.py --sample ewkDM_ttZ_ll_GS_DC2A_0p200000_DC2V_0p200000 #SPLIT20
#python genPostProcessing.py --sample ewkDM_ttZ_ll_GS_DC1A_1p000000 #SPLIT20
#python genPostProcessing.py --sample ewkDM_ttZ_ll_GS_DC1A_0p600000_DC1V_m0p240000_DC2A_0p176700_DC2V_0p176700 #SPLIT20
#python genPostProcessing.py --sample ewkDM_ttZ_ll_GS_DC1A_0p600000_DC1V_m0p240000_DC2A_0p176700_DC2V_m0p176700 #SPLIT20
#python genPostProcessing.py --sample ewkDM_ttZ_ll_GS_DC1A_0p600000_DC1V_m0p240000_DC2A_m0p176700_DC2V_0p176700 #SPLIT20
#python genPostProcessing.py --sample ewkDM_ttZ_ll_GS_DC1A_0p600000_DC1V_m0p240000_DC2A_m0p176700_DC2V_m0p176700 #SPLIT20
#python genPostProcessing.py --sample HEL_UFO_ttZ_ll_GS_cuW_0p100000 #SPLIT20
#python genPostProcessing.py --sample HEL_UFO_ttZ_ll_GS_cuW_0p200000 #SPLIT20
#python genPostProcessing.py --sample HEL_UFO_ttZ_ll_GS_cuW_0p300000 #SPLIT20
#python genPostProcessing.py --sample HEL_UFO_ttZ_ll_GS_cuW_m0p100000 #SPLIT20
#python genPostProcessing.py --sample HEL_UFO_ttZ_ll_GS_cuW_m0p200000 #SPLIT20
#python genPostProcessing.py --sample HEL_UFO_ttZ_ll_GS_cuW_m0p300000 #SPLIT20

#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_1p00_ctZI_1p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m2p00_ctZI_m2p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m2p00_ctZI_m0p80
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m2p00_ctZI_1p60
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_0p00_ctZI_0p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m2p00_ctZI_0p40
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m1p60_ctZI_m1p60
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m0p40_ctZI_2p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m0p80_ctZI_m2p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m0p40_ctZI_m0p40

#python genPostProcessing.py --sample ewkDMGZ_ttgamma_GS #SPLIT10
#python genPostProcessing.py --sample ewkDMGZ_ttgamma_GS_DAG_0p176700_DVG_0p176700 #SPLIT10
#python genPostProcessing.py --sample ewkDMGZ_ttgamma_GS_DAG_0p176700_DVG_m0p176700 #SPLIT10
#python genPostProcessing.py --sample ewkDMGZ_ttgamma_GS_DAG_0p250000 #SPLIT10
#python genPostProcessing.py --sample ewkDMGZ_ttgamma_GS_DAG_0p500000 #SPLIT10
#python genPostProcessing.py --sample ewkDMGZ_ttgamma_GS_DAG_m0p176700_DVG_0p176700 #SPLIT10
#python genPostProcessing.py --sample ewkDMGZ_ttgamma_GS_DAG_m0p176700_DVG_m0p176700 #SPLIT10
#python genPostProcessing.py --sample ewkDMGZ_ttgamma_GS_DAG_m0p250000 #SPLIT10
#python genPostProcessing.py --sample ewkDMGZ_ttgamma_GS_DAG_m0p500000 #SPLIT10
#python genPostProcessing.py --sample ewkDMGZ_ttgamma_GS_DVG_0p250000 #SPLIT10
#python genPostProcessing.py --sample ewkDMGZ_ttgamma_GS_DVG_0p500000 #SPLIT10
#python genPostProcessing.py --sample ewkDMGZ_ttgamma_GS_DVG_m0p250000 #SPLIT10
#python genPostProcessing.py --sample ewkDMGZ_ttgamma_GS_DVG_m0p500000 #SPLIT10

#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m0p80_ctZI_0p40
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m0p80_ctZI_m0p80
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m0p80_ctZI_1p60
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m0p40_ctZI_0p80
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m1p60_ctZI_m0p40
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m0p40_ctZI_m1p60
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_0p00_ctZI_m1p20
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m1p60_ctZI_0p80
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m1p60_ctZI_2p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m1p20_ctZI_m1p20
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m1p20_ctZI_1p20
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_2p00_ctZI_0p80
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_2p00_ctZI_m1p60
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m1p20_ctZI_0p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_2p00_ctZI_m0p40
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_1p60_ctZI_m2p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_1p60_ctZI_1p60
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_0p80_ctZI_2p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_0p40_ctZI_m2p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_1p60_ctZI_0p40
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_0p80_ctZI_m1p60
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_1p60_ctZI_m0p80
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m2p00_ctZI_m1p60
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_0p80_ctZI_0p80
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_1p20_ctZI_m1p20
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_1p20_ctZI_1p20
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_0p40_ctZI_1p60
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_1p20_ctZI_0p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_0p00_ctZI_1p20
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_0p40_ctZI_m0p80
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_0p80_ctZI_m0p40
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_0p40_ctZI_0p40
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_0p00_ctZI_m2p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m0p80_ctZI_m0p40
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m0p80_ctZI_2p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m0p80_ctZI_m1p60
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m0p40_ctZI_0p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m0p40_ctZI_m1p20
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m0p80_ctZI_0p80
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m0p40_ctZI_1p20
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_0p00_ctZI_m0p80
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m2p00_ctZI_2p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m2p00_ctZI_0p80
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m2p00_ctZI_m0p40
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_0p00_ctZI_0p40
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m1p60_ctZI_m1p20
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m1p60_ctZI_0p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m1p60_ctZI_1p20
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_2p00_ctZI_m1p20
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_2p00_ctZI_1p20
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m1p20_ctZI_m2p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m0p80_ctZI_0p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_0p00_ctZI_m1p60
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m1p20_ctZI_m0p80
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m0p40_ctZI_0p40
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m0p40_ctZI_m2p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m0p80_ctZI_m1p20
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m0p40_ctZI_m0p80
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m0p80_ctZI_1p20
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_0p00_ctZI_m0p40
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m0p40_ctZI_1p60
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_2p00_ctZI_0p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m1p20_ctZI_1p60
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m2p00_ctZI_m1p20
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_1p60_ctZI_2p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_1p60_ctZI_m1p60
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m1p20_ctZI_0p40
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_1p20_ctZI_m2p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_1p60_ctZI_0p80
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_1p60_ctZI_m0p40
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_0p40_ctZI_m1p60
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_0p40_ctZI_2p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_1p20_ctZI_m0p80
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_1p20_ctZI_1p60
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_0p80_ctZI_m1p20
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_0p80_ctZI_1p20
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_1p20_ctZI_0p40
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m2p00_ctZI_1p20
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_0p00_ctZI_1p60
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m1p60_ctZI_m2p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_2p00_ctZI_2p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m2p00_ctZI_0p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_0p40_ctZI_m0p40
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_0p40_ctZI_0p80
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_0p80_ctZI_0p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m1p60_ctZI_1p60
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_2p00_ctZI_m0p80
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m1p60_ctZI_m0p80
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m1p60_ctZI_0p40
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_2p00_ctZI_1p60
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m1p20_ctZI_m1p60
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_2p00_ctZI_0p40
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_0p00_ctZI_0p80
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m1p20_ctZI_2p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_2p00_ctZI_m2p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m1p20_ctZI_m0p40
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_1p60_ctZI_m1p20
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_1p20_ctZI_m1p60
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_m1p20_ctZI_0p80
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_1p60_ctZI_1p20
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_1p60_ctZI_0p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_0p80_ctZI_m2p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_1p20_ctZI_2p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_0p80_ctZI_1p60
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_0p80_ctZI_m0p80
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_0p00_ctZI_2p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_1p20_ctZI_m0p40
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_1p20_ctZI_0p80
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_0p40_ctZI_m1p20
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_0p40_ctZI_0p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_0p40_ctZI_1p20
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_ctZ_0p80_ctZI_0p40

#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_m8p00_cpt_m24p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_m8p00_cpt_m14p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_m8p00_cpt_m3p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_m8p00_cpt_7p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_8p00_cpt_m24p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_4p00_cpt_m24p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_4p00_cpt_m14p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_m4p00_cpt_m24p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_4p00_cpt_7p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_8p00_cpt_m14p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_8p00_cpt_7p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_8p00_cpt_m3p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_4p00_cpt_m3p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_0p00_cpt_m24p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_m4p00_cpt_m14p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_m4p00_cpt_m3p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_0p00_cpt_7p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_m4p00_cpt_7p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_0p00_cpt_m14p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_0p00_cpt_m3p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_12p00_cpt_m24p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_16p00_cpt_m24p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_20p00_cpt_m24p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_28p00_cpt_m24p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_20p00_cpt_m14p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_24p00_cpt_m24p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_12p00_cpt_m14p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_m8p00_cpt_m21p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_8p00_cpt_m21p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_12p00_cpt_7p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_16p00_cpt_m3p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_16p00_cpt_m14p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_28p00_cpt_m3p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_4p00_cpt_m21p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_16p00_cpt_7p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_12p00_cpt_m3p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_28p00_cpt_m14p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_20p00_cpt_m3p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_24p00_cpt_m14p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_24p00_cpt_m3p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_24p00_cpt_7p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_28p00_cpt_7p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_20p00_cpt_7p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_4p00_cpt_m10p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_4p00_cpt_10p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_8p00_cpt_10p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_8p00_cpt_0p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_8p00_cpt_m10p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_4p00_cpt_0p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_m8p00_cpt_0p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_m8p00_cpt_m10p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_m8p00_cpt_10p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_m4p00_cpt_m21p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_8p00_cpt_m17p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_4p00_cpt_m17p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_m4p00_cpt_10p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_0p00_cpt_m21p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_4p00_cpt_14p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_4p00_cpt_m7p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_8p00_cpt_3p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_8p00_cpt_14p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_0p00_cpt_m10p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_12p00_cpt_m21p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_16p00_cpt_m21p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_0p00_cpt_10p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_4p00_cpt_3p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_m4p00_cpt_0p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_m4p00_cpt_m10p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_0p00_cpt_0p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_8p00_cpt_m7p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_20p00_cpt_m21p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_28p00_cpt_m21p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_24p00_cpt_m21p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_m8p00_cpt_m17p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_12p00_cpt_m10p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_28p00_cpt_0p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_20p00_cpt_m10p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_24p00_cpt_m10p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_24p00_cpt_0p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_12p00_cpt_10p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_28p00_cpt_m10p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_16p00_cpt_10p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_20p00_cpt_10p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_16p00_cpt_m10p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_28p00_cpt_10p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_20p00_cpt_0p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_24p00_cpt_10p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_m8p00_cpt_3p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_16p00_cpt_0p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_m8p00_cpt_14p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_12p00_cpt_0p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_0p00_cpt_3p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_m8p00_cpt_m7p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_m4p00_cpt_14p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_m4p00_cpt_m17p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_0p00_cpt_m17p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_0p00_cpt_14p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_16p00_cpt_m17p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_12p00_cpt_m17p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_0p00_cpt_m7p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_20p00_cpt_m17p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_m4p00_cpt_m7p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_28p00_cpt_m17p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_24p00_cpt_m17p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_m4p00_cpt_3p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_28p00_cpt_3p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_20p00_cpt_m7p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_24p00_cpt_m7p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_28p00_cpt_m7p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_12p00_cpt_m7p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_12p00_cpt_14p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_24p00_cpt_3p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_28p00_cpt_14p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_16p00_cpt_14p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_20p00_cpt_14p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_16p00_cpt_m7p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_20p00_cpt_3p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_12p00_cpt_3p50
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_24p00_cpt_14p00
#python genPostProcessing.py --sample dim6top_LO_ttZ_ll_cpQM_16p00_cpt_3p50

#python genPostProcessing.py --sample ewkDM_ttZ_ll
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p000000_DC1V_m1p000000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p100000_DC1V_0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p100000_DC1V_0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p100000_DC1V_0p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p100000_DC1V_0p400000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p100000_DC1V_0p500000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p100000_DC1V_0p600000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p100000_DC1V_0p700000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p100000_DC1V_0p800000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p100000_DC1V_0p900000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p100000_DC1V_m0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p100000_DC1V_m0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p100000_DC1V_m0p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p100000_DC1V_m0p400000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p100000_DC1V_m0p500000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p100000_DC1V_m0p600000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p100000_DC1V_m0p700000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p100000_DC1V_m0p800000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p100000_DC1V_m0p900000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p100000_DC1V_m1p000000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p100000_DC1V_m1p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p100000_DC1V_m1p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p100000_DC1V_m1p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p200000_DC1V_0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p200000_DC1V_0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p200000_DC1V_0p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p200000_DC1V_0p400000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p200000_DC1V_0p500000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p200000_DC1V_0p600000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p200000_DC1V_0p700000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p200000_DC1V_0p800000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p200000_DC1V_0p900000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p200000_DC1V_m0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p200000_DC1V_m0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p200000_DC1V_m0p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p200000_DC1V_m0p400000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p200000_DC1V_m0p500000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p200000_DC1V_m0p600000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p200000_DC1V_m0p700000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p200000_DC1V_m0p800000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p200000_DC1V_m0p900000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p200000_DC1V_m1p000000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p200000_DC1V_m1p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p200000_DC1V_m1p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p200000_DC1V_m1p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p300000_DC1V_0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p300000_DC1V_0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p300000_DC1V_0p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p300000_DC1V_0p400000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p300000_DC1V_0p500000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p300000_DC1V_0p600000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p300000_DC1V_0p700000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p300000_DC1V_0p800000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p300000_DC1V_0p900000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p300000_DC1V_m0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p300000_DC1V_m0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p300000_DC1V_m0p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p300000_DC1V_m0p400000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p300000_DC1V_m0p500000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p300000_DC1V_m0p600000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p300000_DC1V_m0p700000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p300000_DC1V_m0p800000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p300000_DC1V_m0p900000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p300000_DC1V_m1p000000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p300000_DC1V_m1p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p300000_DC1V_m1p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p300000_DC1V_m1p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p400000_DC1V_0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p400000_DC1V_0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p400000_DC1V_0p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p400000_DC1V_0p400000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p400000_DC1V_0p500000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p400000_DC1V_0p600000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p400000_DC1V_0p700000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p400000_DC1V_0p800000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p400000_DC1V_0p900000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p400000_DC1V_m0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p400000_DC1V_m0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p400000_DC1V_m0p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p400000_DC1V_m0p400000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p400000_DC1V_m0p500000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p400000_DC1V_m0p600000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p400000_DC1V_m0p700000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p400000_DC1V_m0p800000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p400000_DC1V_m0p900000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p400000_DC1V_m1p000000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p400000_DC1V_m1p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p400000_DC1V_m1p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p400000_DC1V_m1p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p400000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p500000_DC1V_0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p500000_DC1V_0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p500000_DC1V_0p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p500000_DC1V_0p400000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p500000_DC1V_0p500000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p500000_DC1V_0p600000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p500000_DC1V_0p700000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p500000_DC1V_0p800000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p500000_DC1V_0p900000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p500000_DC1V_m0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p500000_DC1V_m0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p500000_DC1V_m0p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p500000_DC1V_m0p400000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p500000_DC1V_m0p500000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p500000_DC1V_m0p600000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p500000_DC1V_m0p700000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p500000_DC1V_m0p800000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p500000_DC1V_m0p900000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p500000_DC1V_m1p000000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p500000_DC1V_m1p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p500000_DC1V_m1p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p500000_DC1V_m1p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p500000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p600000_DC1V_0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p600000_DC1V_0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p600000_DC1V_0p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p600000_DC1V_0p400000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p600000_DC1V_0p500000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p600000_DC1V_0p600000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p600000_DC1V_0p700000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p600000_DC1V_0p800000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p600000_DC1V_0p900000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p600000_DC1V_m0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p600000_DC1V_m0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p600000_DC1V_m0p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p600000_DC1V_m0p400000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p600000_DC1V_m0p500000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p600000_DC1V_m0p600000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p600000_DC1V_m0p700000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p600000_DC1V_m0p800000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p600000_DC1V_m0p900000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p600000_DC1V_m1p000000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p600000_DC1V_m1p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p600000_DC1V_m1p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p600000_DC1V_m1p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p600000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p700000_DC1V_0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p700000_DC1V_0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p700000_DC1V_0p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p700000_DC1V_0p400000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p700000_DC1V_0p500000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p700000_DC1V_0p600000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p700000_DC1V_0p700000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p700000_DC1V_0p800000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p700000_DC1V_0p900000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p700000_DC1V_m0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p700000_DC1V_m0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p700000_DC1V_m0p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p700000_DC1V_m0p400000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p700000_DC1V_m0p500000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p700000_DC1V_m0p600000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p700000_DC1V_m0p700000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p700000_DC1V_m0p800000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p700000_DC1V_m0p900000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p700000_DC1V_m1p000000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p700000_DC1V_m1p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p700000_DC1V_m1p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p700000_DC1V_m1p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p700000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p800000_DC1V_0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p800000_DC1V_0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p800000_DC1V_0p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p800000_DC1V_0p400000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p800000_DC1V_0p500000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p800000_DC1V_0p600000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p800000_DC1V_0p700000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p800000_DC1V_0p800000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p800000_DC1V_0p900000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p800000_DC1V_m0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p800000_DC1V_m0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p800000_DC1V_m0p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p800000_DC1V_m0p400000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p800000_DC1V_m0p500000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p800000_DC1V_m0p600000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p800000_DC1V_m0p700000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p800000_DC1V_m0p800000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p800000_DC1V_m0p900000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p800000_DC1V_m1p000000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p800000_DC1V_m1p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p800000_DC1V_m1p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p800000_DC1V_m1p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p800000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p900000_DC1V_0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p900000_DC1V_0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p900000_DC1V_0p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p900000_DC1V_0p400000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p900000_DC1V_0p500000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p900000_DC1V_0p600000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p900000_DC1V_0p700000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p900000_DC1V_0p800000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p900000_DC1V_0p900000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p900000_DC1V_m0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p900000_DC1V_m0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p900000_DC1V_m0p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p900000_DC1V_m0p400000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p900000_DC1V_m0p500000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p900000_DC1V_m0p600000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p900000_DC1V_m0p700000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p900000_DC1V_m0p800000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p900000_DC1V_m0p900000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p900000_DC1V_m1p000000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p900000_DC1V_m1p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p900000_DC1V_m1p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p900000_DC1V_m1p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_0p900000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p000000_DC1V_0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p000000_DC1V_0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p000000_DC1V_0p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p000000_DC1V_0p400000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p000000_DC1V_0p500000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p000000_DC1V_0p600000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p000000_DC1V_0p700000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p000000_DC1V_0p800000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p000000_DC1V_0p900000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p000000_DC1V_m0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p000000_DC1V_m0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p000000_DC1V_m0p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p000000_DC1V_m0p400000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p000000_DC1V_m0p500000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p000000_DC1V_m0p600000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p000000_DC1V_m0p700000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p000000_DC1V_m0p800000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p000000_DC1V_m0p900000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p000000_DC1V_m1p000000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p000000_DC1V_m1p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p000000_DC1V_m1p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p000000_DC1V_m1p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p000000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p100000_DC1V_0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p100000_DC1V_0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p100000_DC1V_0p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p100000_DC1V_0p400000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p100000_DC1V_0p500000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p100000_DC1V_0p600000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p100000_DC1V_0p700000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p100000_DC1V_0p800000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p100000_DC1V_0p900000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p100000_DC1V_m0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p100000_DC1V_m0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p100000_DC1V_m0p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p100000_DC1V_m0p400000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p100000_DC1V_m0p500000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p100000_DC1V_m0p600000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p100000_DC1V_m0p700000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p100000_DC1V_m0p800000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p100000_DC1V_m0p900000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p100000_DC1V_m1p000000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p100000_DC1V_m1p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p100000_DC1V_m1p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p100000_DC1V_m1p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p200000_DC1V_0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p200000_DC1V_0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p200000_DC1V_0p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p200000_DC1V_0p400000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p200000_DC1V_0p500000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p200000_DC1V_0p600000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p200000_DC1V_0p700000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p200000_DC1V_0p800000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p200000_DC1V_0p900000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p200000_DC1V_m0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p200000_DC1V_m0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p200000_DC1V_m0p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p200000_DC1V_m0p400000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p200000_DC1V_m0p500000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p200000_DC1V_m0p600000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p200000_DC1V_m0p700000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p200000_DC1V_m0p800000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p200000_DC1V_m0p900000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p200000_DC1V_m1p000000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p200000_DC1V_m1p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p200000_DC1V_m1p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p200000_DC1V_m1p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p300000_DC1V_0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p300000_DC1V_0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p300000_DC1V_0p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p300000_DC1V_0p400000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p300000_DC1V_0p500000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p300000_DC1V_0p600000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p300000_DC1V_0p700000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p300000_DC1V_0p800000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p300000_DC1V_0p900000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p300000_DC1V_m0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p300000_DC1V_m0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p300000_DC1V_m0p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p300000_DC1V_m0p400000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p300000_DC1V_m0p500000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p300000_DC1V_m0p600000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p300000_DC1V_m0p700000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p300000_DC1V_m0p800000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p300000_DC1V_m0p900000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p300000_DC1V_m1p000000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p300000_DC1V_m1p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p300000_DC1V_m1p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p300000_DC1V_m1p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p400000_DC1V_0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p400000_DC1V_0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p400000_DC1V_0p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p400000_DC1V_0p400000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p400000_DC1V_0p500000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p400000_DC1V_0p600000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p400000_DC1V_0p700000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p400000_DC1V_0p800000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p400000_DC1V_0p900000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p400000_DC1V_m0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p400000_DC1V_m0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p400000_DC1V_m0p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p400000_DC1V_m0p400000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p400000_DC1V_m0p500000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p400000_DC1V_m0p600000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p400000_DC1V_m0p700000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p400000_DC1V_m0p800000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p400000_DC1V_m0p900000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p400000_DC1V_m1p000000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p400000_DC1V_m1p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p400000_DC1V_m1p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p400000_DC1V_m1p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p400000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p500000_DC1V_0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p500000_DC1V_0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p500000_DC1V_0p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p500000_DC1V_0p400000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p500000_DC1V_0p500000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p500000_DC1V_0p600000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p500000_DC1V_0p700000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p500000_DC1V_0p800000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p500000_DC1V_0p900000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p500000_DC1V_m0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p500000_DC1V_m0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p500000_DC1V_m0p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p500000_DC1V_m0p400000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p500000_DC1V_m0p500000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p500000_DC1V_m0p600000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p500000_DC1V_m0p700000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p500000_DC1V_m0p800000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p500000_DC1V_m0p900000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p500000_DC1V_m1p000000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p500000_DC1V_m1p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p500000_DC1V_m1p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p500000_DC1V_m1p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_1p500000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p100000_DC1V_0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p100000_DC1V_0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p100000_DC1V_0p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p100000_DC1V_0p400000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p100000_DC1V_0p500000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p100000_DC1V_0p600000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p100000_DC1V_0p700000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p100000_DC1V_0p800000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p100000_DC1V_0p900000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p100000_DC1V_m0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p100000_DC1V_m0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p100000_DC1V_m0p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p100000_DC1V_m0p400000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p100000_DC1V_m0p500000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p100000_DC1V_m0p600000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p100000_DC1V_m0p700000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p100000_DC1V_m0p800000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p100000_DC1V_m0p900000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p100000_DC1V_m1p000000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p100000_DC1V_m1p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p100000_DC1V_m1p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p100000_DC1V_m1p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p200000_DC1V_0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p200000_DC1V_0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p200000_DC1V_0p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p200000_DC1V_0p400000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p200000_DC1V_0p500000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p200000_DC1V_0p600000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p200000_DC1V_0p700000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p200000_DC1V_0p800000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p200000_DC1V_0p900000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p200000_DC1V_m0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p200000_DC1V_m0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p200000_DC1V_m0p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p200000_DC1V_m0p400000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p200000_DC1V_m0p500000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p200000_DC1V_m0p600000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p200000_DC1V_m0p700000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p200000_DC1V_m0p800000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p200000_DC1V_m0p900000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p200000_DC1V_m1p000000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p200000_DC1V_m1p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p200000_DC1V_m1p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p200000_DC1V_m1p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p300000_DC1V_0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p300000_DC1V_0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p300000_DC1V_0p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p300000_DC1V_0p400000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p300000_DC1V_0p500000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p300000_DC1V_0p600000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p300000_DC1V_0p700000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p300000_DC1V_0p800000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p300000_DC1V_0p900000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p300000_DC1V_m0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p300000_DC1V_m0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p300000_DC1V_m0p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p300000_DC1V_m0p400000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p300000_DC1V_m0p500000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p300000_DC1V_m0p600000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p300000_DC1V_m0p700000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p300000_DC1V_m0p800000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p300000_DC1V_m0p900000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p300000_DC1V_m1p000000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p300000_DC1V_m1p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p300000_DC1V_m1p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p300000_DC1V_m1p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1A_m0p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1V_0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1V_0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1V_0p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1V_0p400000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1V_0p500000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1V_0p600000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1V_0p700000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1V_0p800000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1V_0p900000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1V_m0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1V_m0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1V_m0p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1V_m0p400000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1V_m0p500000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1V_m0p600000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1V_m0p700000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1V_m0p800000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1V_m0p900000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1V_m1p000000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1V_m1p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1V_m1p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC1V_m1p300000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p050000_DC2V_0p050000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p050000_DC2V_0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p050000_DC2V_0p150000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p050000_DC2V_0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p050000_DC2V_0p250000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p050000_DC2V_m0p050000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p050000_DC2V_m0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p050000_DC2V_m0p150000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p050000_DC2V_m0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p050000_DC2V_m0p250000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p050000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p100000_DC2V_0p050000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p100000_DC2V_0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p100000_DC2V_0p150000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p100000_DC2V_0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p100000_DC2V_0p250000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p100000_DC2V_m0p050000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p100000_DC2V_m0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p100000_DC2V_m0p150000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p100000_DC2V_m0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p100000_DC2V_m0p250000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p150000_DC2V_0p050000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p150000_DC2V_0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p150000_DC2V_0p150000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p150000_DC2V_0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p150000_DC2V_0p250000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p150000_DC2V_m0p050000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p150000_DC2V_m0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p150000_DC2V_m0p150000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p150000_DC2V_m0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p150000_DC2V_m0p250000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p150000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p200000_DC2V_0p050000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p200000_DC2V_0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p200000_DC2V_0p150000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p200000_DC2V_0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p200000_DC2V_0p250000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p200000_DC2V_m0p050000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p200000_DC2V_m0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p200000_DC2V_m0p150000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p200000_DC2V_m0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p200000_DC2V_m0p250000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p250000_DC2V_0p050000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p250000_DC2V_0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p250000_DC2V_0p150000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p250000_DC2V_0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p250000_DC2V_0p250000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p250000_DC2V_m0p050000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p250000_DC2V_m0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p250000_DC2V_m0p150000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p250000_DC2V_m0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p250000_DC2V_m0p250000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_0p250000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p050000_DC2V_0p050000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p050000_DC2V_0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p050000_DC2V_0p150000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p050000_DC2V_0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p050000_DC2V_0p250000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p050000_DC2V_m0p050000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p050000_DC2V_m0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p050000_DC2V_m0p150000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p050000_DC2V_m0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p050000_DC2V_m0p250000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p050000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p100000_DC2V_0p050000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p100000_DC2V_0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p100000_DC2V_0p150000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p100000_DC2V_0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p100000_DC2V_0p250000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p100000_DC2V_m0p050000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p100000_DC2V_m0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p100000_DC2V_m0p150000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p100000_DC2V_m0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p100000_DC2V_m0p250000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p150000_DC2V_0p050000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p150000_DC2V_0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p150000_DC2V_0p150000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p150000_DC2V_0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p150000_DC2V_0p250000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p150000_DC2V_m0p050000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p150000_DC2V_m0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p150000_DC2V_m0p150000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p150000_DC2V_m0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p150000_DC2V_m0p250000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p150000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p200000_DC2V_0p050000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p200000_DC2V_0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p200000_DC2V_0p150000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p200000_DC2V_0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p200000_DC2V_0p250000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p200000_DC2V_m0p050000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p200000_DC2V_m0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p200000_DC2V_m0p150000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p200000_DC2V_m0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p200000_DC2V_m0p250000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p250000_DC2V_0p050000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p250000_DC2V_0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p250000_DC2V_0p150000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p250000_DC2V_0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p250000_DC2V_0p250000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p250000_DC2V_m0p050000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p250000_DC2V_m0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p250000_DC2V_m0p150000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p250000_DC2V_m0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p250000_DC2V_m0p250000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2A_m0p250000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2V_0p050000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2V_0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2V_0p150000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2V_0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2V_0p250000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2V_m0p050000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2V_m0p100000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2V_m0p150000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2V_m0p200000
#python genPostProcessing.py --sample ewkDM_ttZ_ll_DC2V_m0p250000
#python genPostProcessing.py --sample ewkDM_ttZ_ll
