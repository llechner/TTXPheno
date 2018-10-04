
#overwrite="--overwrite"
overwrite=""
small=""
multiplier="4"
level="reco"
version="phase2_v1"

script="NLLPlot"
zRange="0 40" #for NLLPlots

scale="--scale14TeV"
detector="phase2_CMS"

selection='lepSel3-onZ-njet3p-nbjet_medium1p-Zpt0-leptonIso3'
#nohup krenew -t -K 10 -- bash -c "python ${script}.py ${small} --process ttZ_3l --selection ${selection} --sample fwlite_ttZ_ll_LO_order2_15weights_ref --level ${level} --version ${version} --variables cpQM cpt --zRange ${zRange} --luminosity 3000 --binning 31 -8 28 31 -24 14 --fit BestFit --smooth --contours --cores 6 ${overwrite} --binMultiplier ${multiplier} ${scale} --detector ${detector}" &
nohup krenew -t -K 10 -- bash -c "python ${script}.py ${small} --process ttZ_3l --selection ${selection} --sample fwlite_ttZ_ll_LO_order2_15weights_ref --level ${level} --version ${version}_semizoom --variables cpQM cpt --zRange ${zRange} --luminosity 3000 --binning 31 -8 12 31 -8 12 --fit BestFit --smooth --contours --cores 6 ${overwrite} --binMultiplier ${multiplier} ${scale} --detector ${detector}" &
#nohup krenew -t -K 10 -- bash -c "python ${script}.py ${small} --process ttZ_3l --selection ${selection} --sample fwlite_ttZ_ll_LO_order2_15weights_ref --level ${level} --version ${version}_zoom --variables cpQM cpt --zRange ${zRange} --luminosity 3000 --binning 31 -4 4 31 -4 4 --fit BestFit --smooth --contours --cores 6 ${overwrite} --binMultiplier ${multiplier} ${scale} --detector ${detector}" &
#nohup krenew -t -K 10 -- bash -c "python ${script}.py ${small} --process ttZ_3l --selection ${selection} --sample fwlite_ttZ_ll_LO_order2_15weights_ref --level ${level} --version ${version} --variables ctZ ctZI --zRange ${zRange} --luminosity 3000 --binning 31 -1.5 1.5 31 -1.5 1.5 --fit BestFit --smooth --contours --cores 6 ${overwrite} --binMultiplier ${multiplier} ${scale} --detector ${detector}" &

#nohup krenew -t -K 10 -- bash -c "python ${script}.py ${small} --process ttZ_3l --selection ${selection} --sample fwlite_ttZ_ll_LO_order2_15weights_ref --level ${level} --version ${version} --variables cpQM cpt --zRange ${zRange} --luminosity 3000 --binning 15 -8 28 15 -24 14 --fit BestFit --smooth --contours --cores 2 ${overwrite} --binMultiplier ${multiplier} ${scale} --detector ${detector}" &
#nohup krenew -t -K 10 -- bash -c "python ${script}.py ${small} --process ttZ_3l --selection ${selection} --sample fwlite_ttZ_ll_LO_order2_15weights_ref --level ${level} --version ${version}_zoom --variables cpQM cpt --zRange ${zRange} --luminosity 3000 --binning 15 -4 4 15 -4 4 --fit BestFit --smooth --contours --cores 2 ${overwrite} --binMultiplier ${multiplier} ${scale} --detector ${detector}" &
#nohup krenew -t -K 10 -- bash -c "python ${script}.py ${small} --process ttZ_3l --selection ${selection} --sample fwlite_ttZ_ll_LO_order2_15weights_ref --level ${level} --version ${version} --variables ctZ ctZI --zRange ${zRange} --luminosity 3000 --binning 15 -1.5 1.5 15 -1.5 1.5 --fit BestFit --smooth --contours --cores 2 ${overwrite} --binMultiplier ${multiplier} ${scale} --detector ${detector}" &

