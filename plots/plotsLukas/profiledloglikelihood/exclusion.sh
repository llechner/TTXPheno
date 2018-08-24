small=""
level="reco"
version="test"
script="ExclusionPlot.py"
#script="NLLPlot.py"

selection='lepSel3-onZ-njet3p-nbjet1p-Zpt0-leptonIso3'
nohup krenew -t -K 10 -- bash -c "python ${script} ${small} --process ttZ_3l --selection ${selection} --sample fwlite_ttZ_ll_LO_order2_15weights_ref --level ${level} --version ${version} --variables cpQM cpt --zRange 0 40 --luminosity 150 --binning 11 -8 28 11 -24 14 --fit BestFit --cores 4" &
nohup krenew -t -K 10 -- bash -c "python ${script} ${small} --process ttZ_3l --selection ${selection} --sample fwlite_ttZ_ll_LO_order2_15weights_ref --level ${level} --version ${version} --variables ctZ ctZI --zRange 0 40 --luminosity 150 --binning 11 -1.5 1.5 11 -1.5 1.5 --fit BestFit --cores 4" &

