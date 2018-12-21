
#overwrite="--overwrite"
overwrite=""
small=""
multiplier="10"
level="reco"
#version="paper_250818_3000"
version="paper_250818"

#script="ExclusionPlot"
#zRange="0.1 100" #for Exclusion plots

script="NLLPlot"
zRange="0 40" #for NLLPlots

scale="--scale 3000"
#scale=""

selection='lepSel3-onZ-njet3p-nbjet1p-Zpt0-leptonIso3'
nohup krenew -t -K 10 -- bash -c "python ${script}.py ${small} --process ttZ_3l --selection ${selection} --sample fwlite_ttZ_ll_LO_order2_15weights_ref --level ${level} --version ${version} --variables cpQM cpt --zRange ${zRange} --luminosity 150 --binning 31 -8 28 31 -24 14 --fit BestFit --smooth --contours --cores 2 ${overwrite} --binMultiplier ${multiplier} ${scale}" &
nohup krenew -t -K 10 -- bash -c "python ${script}.py ${small} --process ttZ_3l --selection ${selection} --sample fwlite_ttZ_ll_LO_order2_15weights_ref --level ${level} --version ${version}_zoom --variables cpQM cpt --zRange ${zRange} --luminosity 150 --binning 31 -4 4 31 -4 4 --fit BestFit --smooth --contours --cores 6 ${overwrite} --binMultiplier ${multiplier} ${scale}" &
#nohup krenew -t -K 10 -- bash -c "python ${script}.py ${small} --process ttZ_3l --selection ${selection} --sample fwlite_ttZ_ll_LO_order2_15weights_ref --level ${level} --version ${version} --variables ctZ ctZI --zRange ${zRange} --luminosity 150 --binning 31 -1.5 1.5 31 -1.5 1.5 --fit BestFit --smooth --contours --cores 2 ${overwrite} --binMultiplier ${multiplier} ${scale}" &

#CuB CuW
#nohup krenew -t -K 10 -- bash -c "python ${script}.py ${small} --process ttZ_3l --selection ${selection} --sample fwlite_ttZ_ll_LO_order2_15weights_ref --level ${level} --version ${version} --variables cuB cuW --zRange ${zRange} --luminosity 150 --binning 31 -2 2 31 -2 2 --fit BestFit --smooth --contours --cores 2 ${overwrite} --binMultiplier ${multiplier}" &


selection='lepSel1-gammapt40-njet3p-nbjet1p-relIso0to0.12-met40-leptonIso1'
#nohup krenew -t -K 10 -- bash -c "python ${script}.py ${small} --process ttgamma_1l --selection ${selection} --sample fwlite_ttgammaLarge_LO_order2_15weights_ref --level ${level} --version ${version} --variables ctZ ctZI --zRange ${zRange} --luminosity 150 --binning 31 -0.5 0.5 31 -0.5 0.5 --smooth --contours --cores 2 --fit BestFit ${overwrite} --binMultiplier ${multiplier} ${scale}" &

#CuB CuW
#nohup krenew -t -K 10 -- bash -c "python ${script}.py ${small} --process ttgamma_1l --selection ${selection} --sample fwlite_ttgammaLarge_LO_order2_15weights_ref --level ${level} --version ${version} --variables cuB cuW --zRange ${zRange} --luminosity 150 --binning 31 -2 2 31 -2 2 --smooth --contours --cores 2 --fit BestFit ${overwrite} --binMultiplier ${multiplier}" &


selection='lepSel2-gammapt40-njet2p-nbjet1p-relIso0to0.12-met40-leptonIso2'
#nohup krenew -t -K 10 -- bash -c "python ${script}.py ${small} --process ttgamma_2l --selection ${selection} --sample fwlite_ttgammaLarge_LO_order2_15weights_ref --level ${level} --version ${version} --variables ctZ ctZI --zRange ${zRange} --luminosity 150 --binning 31 -0.5 0.5 31 -0.5 0.5 --smooth --contours --cores 3 --fit BestFit ${overwrite} --binMultiplier ${multiplier} ${scale}" &

#CuB CuW
#nohup krenew -t -K 10 -- bash -c "python ${script}.py ${small} --process ttgamma_2l --selection ${selection} --sample fwlite_ttgammaLarge_LO_order2_15weights_ref --level ${level} --version ${version} --variables cuB cuW --zRange ${zRange} --luminosity 150 --binning 31 -2 2 31 -2 2 --smooth --contours --cores 3 --fit BestFit ${overwrite} --binMultiplier ${multiplier}" &

#combined
#nohup krenew -t -K 10 -- bash -c "python ${script}_combined.py ${small} --level ${level} --version ${version} --variables ctZ ctZI --binning 31 -0.3 0.3 31 -0.3 0.3 --luminosity 150 --fit BestFit --zRange ${zRange} --contours --smooth --cores 6 ${overwrite} --binMultiplier ${multiplier} ${scale}" &
#nohup krenew -t -K 10 -- bash -c "python ${script}_combined.py ${small} --level ${level} --version ${version} --variables ctZ ctZI --binning 11 -0.3 0.3 11 -0.3 0.3 --luminosity 150 --fit BestFit --zRange ${zRange} --contours --smooth --cores 2 ${overwrite} --binMultiplier ${multiplier} --small" &

#CuB CuW
#nohup krenew -t -K 10 -- bash -c "python ${script}_combined.py ${small} --level ${level} --version ${version} --variables cuB cuW --binning 31 -2 2 31 -2 2 --luminosity 150 --fit BestFit --zRange ${zRange} --contours --smooth --cores 6 ${overwrite} --binMultiplier ${multiplier}" &


selection='lepSel3-onZ-njet3p-nbjet1p-Zpt0-leptonIso3'
#nohup krenew -t -K 10 -- bash -c "python ${script}.py ${small} --process ttZ_3l --selection ${selection} --sample fwlite_ttZ_ll_LO_order2_15weights_ref --level ${level} --version ${version} --variables ctZ ctZI --zRange ${zRange} --luminosity 3000 --binning 31 -1.5 1.5 31 -1.5 1.5 --fit BestFit --smooth --contours --cores 6 ${overwrite} --binMultiplier ${multiplier} ${scale}" &

selection='lepSel1-gammapt40-njet3p-nbjet1p-relIso0to0.12-met40-leptonIso1'
#nohup krenew -t -K 10 -- bash -c "python ${script}.py ${small} --process ttgamma_1l --selection ${selection} --sample fwlite_ttgammaLarge_LO_order2_15weights_ref --level ${level} --version ${version} --variables ctZ ctZI --zRange ${zRange} --luminosity 3000 --binning 31 -0.5 0.5 31 -0.5 0.5 --smooth --contours --cores 6 --fit BestFit ${overwrite} --binMultiplier ${multiplier} ${scale}" &

selection='lepSel2-gammapt40-njet2p-nbjet1p-relIso0to0.12-met40-leptonIso2'
#nohup krenew -t -K 10 -- bash -c "python ${script}.py ${small} --process ttgamma_2l --selection ${selection} --sample fwlite_ttgammaLarge_LO_order2_15weights_ref --level ${level} --version ${version} --variables ctZ ctZI --zRange ${zRange} --luminosity 3000 --binning 31 -0.5 0.5 31 -0.5 0.5 --smooth --contours --cores 6 --fit BestFit ${overwrite} --binMultiplier ${multiplier} ${scale}" &

#nohup krenew -t -K 10 -- bash -c "python ${script}_combined.py ${small} --level ${level} --version ${version} --variables ctZ ctZI --binning 31 -0.3 0.3 31 -0.3 0.3 --luminosity 3000 --fit BestFit --zRange ${zRange} --contours --smooth --cores 6 ${overwrite} --binMultiplier ${multiplier} ${scale}" &
