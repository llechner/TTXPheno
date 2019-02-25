version="note_v23_ARC"
detector='phase2_CMS'

nohup krenew -t -K 10 -- bash -c "python skim_plots.py --processFile ttZ_3l_paper --luminosity 3000 --version ${version} --level reco --sample fwlite_ttZ_ll_LO_order2_15weights_ref --order 2 --selection lepSel3-onZ-njet3p-nbjet1p-Zpt0-leptonIso3 --parameters ctZ 2 ctZI 2 --backgrounds --detector ${detector} --noninfoSignal --scaleLumi --scale14TeV" &
nohup krenew -t -K 10 -- bash -c "python skim_plots.py --processFile ttZ_3l_paper --luminosity 3000 --version ${version} --level reco --sample fwlite_ttZ_ll_LO_order2_15weights_ref --order 2 --selection lepSel3-onZ-njet3p-nbjet1p-Zpt200-leptonIso3 --parameters ctZ 2 ctZI 2 --backgrounds --detector ${detector} --noninfoSignal --scaleLumi --scale14TeV --wideLegendStyle" &

version="tmp"
nohup krenew -t -K 10 -- bash -c "python skim_plots.py --processFile ttZ_3l_paper --luminosity 3000 --version ${version} --level reco --sample fwlite_ttZ_ll_LO_order2_15weights_ref --order 2 --selection lepSel3-onZ-njet3p-nbjet1p-Zpt0-leptonIso3 --parameters ctZ 2 ctZI 2 --backgrounds --detector ${detector} --noninfoSignal --scaleLumi --scale14TeV --small" &
nohup krenew -t -K 10 -- bash -c "python skim_plots.py --processFile ttZ_3l_paper --luminosity 3000 --version ${version} --level reco --sample fwlite_ttZ_ll_LO_order2_15weights_ref --order 2 --selection lepSel3-onZ-njet3p-nbjet1p-Zpt200-leptonIso3 --parameters ctZ 2 ctZI 2 --backgrounds --detector ${detector} --noninfoSignal --scaleLumi --scale14TeV --wideLegendStyle --small" &

#nohup krenew -t -K 10 -- bash -c "python skim_plots.py --processFile ttZ_3l_paper --luminosity 3000 --version ${version}_semizoom --level reco --sample fwlite_ttZ_ll_LO_order2_15weights_ref --order 2 --selection lepSel3-onZ-njet3p-nbjet1p-Zpt0-leptonIso3 --parameters ctZ 2.4 ctZI 2.4 --backgrounds --detector ${detector} --noninfoSignal --scaleLumi --scale14TeV --wideLegendStyle" &


#nohup krenew -t -K 10 -- bash -c "python skim_plots.py --processFile ttZ_3l_gen --luminosity 3000 --version ${version} --level gen --sample fwlite_ttZ_ll_LO_order2_15weights_ref --order 2 --selection lepSel3-onZ-njet3p-nbjet1p-Zpt0-leptonIso3 --parameters ctZ 2.4 ctZI 2.4 --backgrounds --detector ${detector} --noninfoSignal --scaleLumi --scale14TeV --wideLegendStyle" &

#echo "python skim_plots.py --processFile ttZ_3l_paper --luminosity 3000 --version ${version} --level reco --sample fwlite_ttZ_ll_LO_order2_15weights_ref --order 2 --selection lepSel3-onZ-njet3p-nbjet1p-Zpt0-leptonIso3 --parameters ctZ 2.4 ctZI 2.4 --backgrounds --detector ${detector} --noninfoSignal --scaleLumi --scale14TeV" &

#nohup krenew -t -K 10 -- bash -c "python skim_plots.py --processFile ttZ_3l_paper --luminosity 3000 --version ${version} --level reco --sample fwlite_ttZ_ll_LO_order2_15weights_ref --order 2 --selection lepSel3-onZ-njet3p-nbjet1p-Zpt300-leptonIso3 --parameters ctZ 2.4 ctZI 2.4 --backgrounds --detector ${detector} --noninfoSignal --scaleLumi --scale14TeV" &
#nohup krenew -t -K 10 -- bash -c "python skim_plots.py --processFile ttZ_3l_paper --luminosity 3000 --version ${version} --level reco --sample fwlite_ttZ_ll_LO_order2_15weights_ref --order 2 --selection lepSel3-onZ-njet3p-nbjet1p-Zpt300-leptonIso3 --parameters ctZ 2.4 ctZI 2.4 --backgrounds --detector ${detector} --noninfoSignal --scaleLumi --scale14TeV --small" &
