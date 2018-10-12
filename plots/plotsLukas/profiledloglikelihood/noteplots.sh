
overwrite="--overwrite"
#overwrite=""
small=""
multiplier="5"
level="reco"
version="note_v7"

script="NLLPlot"
zRange="0 40" #for NLLPlots

scale="--scale14TeV"
detector="phase2_CMS"

nonPrompt="--addNonPrompt"
#nonPrompt=""

#removeCards="--removeCardFiles"
removeCards=""

bestFit="--bestFit"
#bestFit=""

combine="--useCombine"
#combine=""

#fitOnly="--fitOnly"
fitOnly=""

selection='lepSel3-onZ-njet3p-nbjet1p-Zpt0-leptonIso3'
nohup krenew -t -K 10 -- bash -c "python ${script}.py ${small} --process ttZ_3l --selection ${selection} --sample fwlite_ttZ_ll_LO_order2_15weights_ref --level ${level} --version ${version} --variables cpQM cpt --zRange ${zRange} --luminosity 3000 --binning 31 -8 28 31 -24 14 --smooth --contours --cores 6 ${overwrite} --binMultiplier ${multiplier} ${scale} --detector ${detector} ${nonPrompt} ${combine} ${removeCards} ${bestFit} ${fitOnly}" &
#nohup krenew -t -K 10 -- bash -c "python ${script}.py ${small} --process ttZ_3l --selection ${selection} --sample fwlite_ttZ_ll_LO_order2_15weights_ref --level ${level} --version ${version}_semizoom --variables cpQM cpt --zRange ${zRange} --luminosity 3000 --binning 31 -3 7 31 -3 7 --smooth --contours --cores 6 ${overwrite} --binMultiplier ${multiplier} ${scale} --detector ${detector} ${nonPrompt} ${combine} ${removeCards} ${bestFit} ${fitOnly}" &
#nohup krenew -t -K 10 -- bash -c "python ${script}.py ${small} --process ttZ_3l --selection ${selection} --sample fwlite_ttZ_ll_LO_order2_15weights_ref --level ${level} --version ${version}_zoom --variables cpQM cpt --zRange ${zRange} --luminosity 3000 --binning 31 -2 2 31 -2 2 --smooth --contours --cores 6 ${overwrite} --binMultiplier ${multiplier} ${scale} --detector ${detector} ${nonPrompt} ${combine} ${removeCards} ${bestFit} ${fitOnly}" &
#nohup krenew -t -K 10 -- bash -c "python ${script}.py ${small} --process ttZ_3l --selection ${selection} --sample fwlite_ttZ_ll_LO_order2_15weights_ref --level ${level} --version ${version} --variables ctZ ctZI --zRange ${zRange} --luminosity 3000 --binning 31 -1.0 1.0 31 -1.0 1.0 --smooth --contours --cores 6 ${overwrite} --binMultiplier ${multiplier} ${scale} --detector ${detector} ${nonPrompt} ${combine} ${removeCards} ${bestFit} ${fitOnly}" &

#nohup krenew -t -K 10 -- bash -c "python ${script}.py ${small} --process ttZ_3l --selection ${selection} --sample fwlite_ttZ_ll_LO_order2_15weights_ref --level ${level} --version ${version} --variables cpQM cpt --zRange ${zRange} --luminosity 3000 --binning 21 -8 28 21 -24 14 --smooth --contours --cores 6 ${overwrite} --binMultiplier ${multiplier} ${scale} --detector ${detector} ${nonPrompt} ${combine} ${removeCards} ${bestFit} ${fitOnly}" &
#nohup krenew -t -K 10 -- bash -c "python ${script}.py ${small} --process ttZ_3l --selection ${selection} --sample fwlite_ttZ_ll_LO_order2_15weights_ref --level ${level} --version ${version}_semizoom --variables cpQM cpt --zRange ${zRange} --luminosity 3000 --binning 21 -8 12 21 -8 12 --smooth --contours --cores 6 ${overwrite} --binMultiplier ${multiplier} ${scale} --detector ${detector} ${nonPrompt} ${combine} ${removeCards} ${bestFit} ${fitOnly}" &
#nohup krenew -t -K 10 -- bash -c "python ${script}.py ${small} --process ttZ_3l --selection ${selection} --sample fwlite_ttZ_ll_LO_order2_15weights_ref --level ${level} --version ${version}_zoom --variables cpQM cpt --zRange ${zRange} --luminosity 3000 --binning 21 -4 4 21 -4 4 --smooth --contours --cores 6 ${overwrite} --binMultiplier ${multiplier} ${scale} --detector ${detector} ${nonPrompt} ${combine} ${removeCards} ${bestFit} ${fitOnly}" &
#nohup krenew -t -K 10 -- bash -c "python ${script}.py ${small} --process ttZ_3l --selection ${selection} --sample fwlite_ttZ_ll_LO_order2_15weights_ref --level ${level} --version ${version} --variables ctZ ctZI --zRange ${zRange} --luminosity 3000 --binning 21 -1.5 1.5 21 -1.5 1.5 --smooth --contours --cores 6 ${overwrite} --binMultiplier ${multiplier} ${scale} --detector ${detector} ${nonPrompt} ${combine} ${removeCards} ${bestFit} ${fitOnly}" &

#echo "python ${script}.py ${small} --process ttZ_3l --selection ${selection} --sample fwlite_ttZ_ll_LO_order2_15weights_ref --level ${level} --version ${version}_semizoom --variables cpQM cpt --zRange ${zRange} --luminosity 3000 --binning 21 -8 12 21 -8 12 --smooth --contours --cores 6 ${overwrite} --binMultiplier ${multiplier} ${scale} --detector ${detector} ${nonPrompt} ${combine} ${removeCards} ${bestFit} ${fitOnly}" &
#echo "python ${script}.py ${small} --process ttZ_3l --selection ${selection} --sample fwlite_ttZ_ll_LO_order2_15weights_ref --level ${level} --version ${version} --variables ctZ ctZI --zRange ${zRange} --luminosity 3000 --binning 21 -1.5 1.5 21 -1.5 1.5 --smooth --contours --cores 6 ${overwrite} --binMultiplier ${multiplier} ${scale} --detector ${detector} ${nonPrompt} ${combine} ${removeCards} ${bestFit} ${fitOnly}" &

#nohup krenew -t -K 10 -- bash -c "python ${script}.py ${small} --process ttZ_3l --selection ${selection} --sample fwlite_ttZ_ll_LO_order2_15weights_ref --level ${level} --version ${version}_semizoom --variables cpQM cpt --zRange ${zRange} --luminosity 3000 --binning 21 -2 6 21 -2 6 --smooth --contours --cores 6 ${overwrite} --binMultiplier ${multiplier} ${scale} --detector ${detector} ${nonPrompt} ${combine} ${removeCards} ${bestFit} ${fitOnly}" &
