version="paper_test"
#version="v18"
detector='CMS'
#detector='ATLAS'

#nohup krenew -t -K 10 -- bash -c "python skim_plots.py --processFile ttZ_3l_paper --luminosity 150 --version ${version} --level reco   --scaleLumi --sample fwlite_ttZ_ll_LO_order2_15weights_ref --order 2 --selection lepSel3-onZ-njet3p-nbjet1p-Zpt0-leptonIso3 --parameters cpQM 7 cpt 8 ctZ 6 ctZI 6 --backgrounds --detector ${detector} --noninfoSignal" &
#nohup krenew -t -K 10 -- bash -c "python skim_plots.py --processFile ttZ_3l_paper --luminosity 150 --version ${version} --level reco   --scaleLumi --sample fwlite_ttZ_ll_LO_order2_15weights_ref --order 2 --selection lepSel3-onZ-njet3p-nbjet1p-Zpt200-leptonIso3 --parameters cpQM 7 cpt 8 ctZ 6 ctZI 6 --backgrounds --detector ${detector} --noninfoSignal" &
#nohup krenew -t -K 10 -- bash -c "python skim_plots.py --processFile ttgamma_1l_paper --luminosity 150 --version ${version} --level reco   --scaleLumi --sample fwlite_ttgammaLarge_LO_order2_15weights_ref --order 2 --selection lepSel1-gammapt40-njet3p-nbjet1p-relIso0to0.12-met40-leptonIso1 --backgrounds --parameters cpQM 7 cpt 8 ctZ 6 ctZI 6 --detector ${detector} --noninfoSignal" &
#nohup krenew -t -K 10 -- bash -c "python skim_plots.py --processFile ttgamma_2l_paper --luminosity 150 --version ${version} --level reco   --scaleLumi --sample fwlite_ttgammaLarge_LO_order2_15weights_ref --order 2 --selection lepSel2-gammapt40-njet2p-nbjet1p-relIso0to0.12-met40-leptonIso2 --backgrounds --parameters cpQM 7 cpt 8 ctZ 6 ctZI 6 --detector ${detector} --noninfoSignal" &

#nohup krenew -t -K 10 -- bash -c "python skim_plots.py --processFile ttZ_3l_paper --luminosity 150 --version ${version} --level reco   --scaleLumi --sample fwlite_ttZ_ll_LO_order2_15weights_ref --order 2 --selection lepSel3-onZ-njet3p-nbjet1p-Zpt0-leptonIso3 --parameters cpQM 7 cpt 8 ctZ 6 ctZI 6 --backgrounds --addFisherInformation --binThreshold 100 --variables cpt --detector ${detector} --noninfoSignal" &
#nohup krenew -t -K 10 -- bash -c "python skim_plots.py --processFile ttZ_3l_paper --luminosity 150 --version ${version} --level reco   --scaleLumi --sample fwlite_ttZ_ll_LO_order2_15weights_ref --order 2 --selection lepSel3-onZ-njet3p-nbjet1p-Zpt0-leptonIso3 --parameters cpQM 7 cpt 8 ctZ 6 ctZI 6 --backgrounds --addFisherInformation --binThreshold 100 --variables cpQM --detector ${detector} --noninfoSignal" &
#nohup krenew -t -K 10 -- bash -c "python skim_plots.py --processFile ttZ_3l_paper --luminosity 150 --version ${version} --level reco   --scaleLumi --sample fwlite_ttZ_ll_LO_order2_15weights_ref --order 2 --selection lepSel3-onZ-njet3p-nbjet1p-Zpt200-leptonIso3 --parameters cpQM 7 cpt 8 ctZ 6 ctZI 6 --backgrounds --addFisherInformation --binThreshold 100 --variables cpt --detector ${detector} --noninfoSignal" &
#nohup krenew -t -K 10 -- bash -c "python skim_plots.py --processFile ttZ_3l_paper --luminosity 150 --version ${version} --level reco   --scaleLumi --sample fwlite_ttZ_ll_LO_order2_15weights_ref --order 2 --selection lepSel3-onZ-njet3p-nbjet1p-Zpt200-leptonIso3 --parameters cpQM 7 cpt 8 ctZ 6 ctZI 6 --backgrounds --addFisherInformation --binThreshold 100 --variables cpQM --detector ${detector} --noninfoSignal" &

#nohup krenew -t -K 10 -- bash -c "python skim_plots.py --processFile ttZ_3l_paper --luminosity 150 --version ${version} --level reco  --sample fwlite_ttZ_ll_LO_order2_15weights_ref --order 2 --selection lepSel3-onZ-njet3p-nbjet1p-Zpt0-leptonIso3 --parameters cpt 8 --addFisherInformation --binThreshold 100 --variables cpt --detector ${detector} --addFisherInformationBackground" &
#nohup krenew -t -K 10 -- bash -c "python skim_plots.py --processFile ttZ_3l_paper --luminosity 150 --version ${version} --level reco  --sample fwlite_ttZ_ll_LO_order2_15weights_ref --order 2 --selection lepSel3-onZ-njet3p-nbjet1p-Zpt0-leptonIso3 --parameters cpQM 7 --addFisherInformation --binThreshold 100 --variables cpQM --detector ${detector} --addFisherInformationBackground" &
#nohup krenew -t -K 10 -- bash -c "python skim_plots.py --processFile ttZ_3l_paper --luminosity 150 --version ${version} --level reco  --sample fwlite_ttZ_ll_LO_order2_15weights_ref --order 2 --selection lepSel3-onZ-njet3p-nbjet1p-Zpt0-leptonIso3 --parameters cpt 8 --addFisherInformation --binThreshold 100 --variables cpt --detector ${detector} --addFisherInformationBackground --scaleLumi" &
#nohup krenew -t -K 10 -- bash -c "python skim_plots.py --processFile ttZ_3l_paper --luminosity 150 --version ${version} --level reco  --sample fwlite_ttZ_ll_LO_order2_15weights_ref --order 2 --selection lepSel3-onZ-njet3p-nbjet1p-Zpt0-leptonIso3 --parameters cpQM 7 --addFisherInformation --binThreshold 100 --variables cpQM --detector ${detector} --addFisherInformationBackground --scaleLumi" &

nohup krenew -t -K 10 -- bash -c "python skim_plots.py --processFile ttZ_3l_paper --luminosity 150 --version ${version} --level reco --sample fwlite_ttZ_ll_LO_order2_15weights_ref --order 2 --selection lepSel3-onZ-njet3p-nbjet1p-Zpt0-leptonIso3 --parameters cpt 8 --backgrounds --detector ${detector} --noninfoSignal" &
nohup krenew -t -K 10 -- bash -c "python skim_plots.py --processFile ttZ_3l_paper --luminosity 150 --version ${version} --level reco --sample fwlite_ttZ_ll_LO_order2_15weights_ref --order 2 --selection lepSel3-onZ-njet3p-nbjet1p-Zpt200-leptonIso3 --parameters cpQM 7 --backgrounds --detector ${detector} --noninfoSignal" &