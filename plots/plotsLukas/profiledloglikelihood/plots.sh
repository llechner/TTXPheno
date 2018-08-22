
#small="--small"
small=""
level="reco"
#version="TTXPheno_08082018"
version="nll"
#script="ExclusionPlot.py"
script="NLLPlot.py"

selection='lepSel3-onZ-njet3p-nbjet1p-Zpt0-mll12'
#python ${script} ${small} --process ttZ_3l --selection ${selection} --sample fwlite_ttZ_ll_LO_order2_15weights_ref --level ${level} --version ${version} --variables cpQM cpt --zRange 0 40 --luminosity 150 --binning 21 -8 28 21 -24 14
#python ${script} ${small} --process ttZ_3l --selection ${selection} --sample fwlite_ttZ_ll_LO_order2_15weights_ref --level ${level} --version ${version} --variables ctZ ctZI --zRange 0 40 --luminosity 150 --binning 21 -2 2 21 -2 2

selection='lepSel1-gammapt40-njet3p-nbjet1p-relIso0to0.12-met40'
nohup krenew -t -K 10 -- bash -c "python ${script} ${small} --process ttgamma_1l --selection ${selection} --sample fwlite_ttgamma_LO_order2_15weights_ref --level ${level} --version ${version} --variables cpQM cpt --zRange 0 40 --luminosity 150 --binning 21 -8 28 21 -24 14 --smooth --contours --cores 2" &
nohup krenew -t -K 10 -- bash -c "python ${script} ${small} --process ttgamma_1l --selection ${selection} --sample fwlite_ttgamma_LO_order2_15weights_ref --level ${level} --version ${version} --variables ctZ ctZI --zRange 0 40 --luminosity 150 --binning 21 -2 2 21 -2 2 --smooth --contours --cores 2" &
nohup krenew -t -K 10 -- bash -c "python ${script} ${small} --process ttgamma_1l --selection ${selection} --sample fwlite_ttgamma_LO_order2_15weights_ref --level ${level} --version ${version} --variables ctW ctWI --zRange 0 40 --luminosity 150 --binning 21 -2 2 21 -2 2 --smooth --contours --cores 2" &

selection='lepSel2-gammapt40-njet2p-nbjet1p-relIso0to0.12-met40'
#nohup krenew -t -K 10 -- bash -c "python ${script} ${small} --process ttgamma_2l --selection ${selection} --sample fwlite_ttgamma_LO_order2_15weights_ref --level ${level} --version ${version} --variables cpQM cpt --zRange 0 40 --luminosity 150 --binning 21 -8 28 21 -24 14 --smooth --contours --cores 2" &
#nohup krenew -t -K 10 -- bash -c "python ${script} ${small} --process ttgamma_2l --selection ${selection} --sample fwlite_ttgamma_LO_order2_15weights_ref --level ${level} --version ${version} --variables ctZ ctZI --zRange 0 40 --luminosity 150 --binning 21 -2 2 21 -2 2 --smooth --contours --cores 2" &
#nohup krenew -t -K 10 -- bash -c "python ${script} ${small} --process ttgamma_2l --selection ${selection} --sample fwlite_ttgamma_LO_order2_15weights_ref --level ${level} --version ${version} --variables ctW ctWI --zRange 0 40 --luminosity 150 --binning 21 -2 2 21 -2 2 --smooth --contours --cores 2" &

