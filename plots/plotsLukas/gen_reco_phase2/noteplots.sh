version="note_v7"
detector='phase2_CMS'

nohup krenew -t -K 10 -- bash -c "python skim_plots.py --processFile ttZ_3l_paper --luminosity 3000 --version ${version} --level reco --sample fwlite_ttZ_ll_LO_order2_15weights_ref --order 2 --selection lepSel3-onZ-njet3p-nbjet1p-Zpt0-leptonIso3 --parameters ctZ 6 ctZI 6 cpt 8 cpQM 7 --backgrounds --detector ${detector} --noninfoSignal --scaleLumi --scale14TeV" &
nohup krenew -t -K 10 -- bash -c "python skim_plots.py --processFile ttZ_3l_paper --luminosity 3000 --version ${version} --level reco --sample fwlite_ttZ_ll_LO_order2_15weights_ref --order 2 --selection lepSel3-onZ-njet3p-nbjet1p-Zpt200-leptonIso3 --parameters ctZ 6 ctZI 6 cpt 8 cpQM 7 --backgrounds --detector ${detector} --noninfoSignal --scaleLumi --scale14TeV" &
