
small=""
level="reco"
version="note_v22"

scale="--scale14TeV"
detector="phase2_CMS"

nonPrompt="--addNonPrompt"
#nonPrompt=""

#combine="--combineTTX"
combine=""

#addOthers="--addOthers"
addOthers=""

#selection='lepSel3-onZ-njet3p-nbjet1p-Zpt0-leptonIso3'
#nohup krenew -t -K 10 -- bash -c "python signalRegion.py ${small} --process ttZ --selection ${selection} --sample fwlite_ttZ_ll_LO_order2_15weights_ref --level ${level} --version ${version} --parameters ctZ 2.4 --luminosity 3000 ${scale} --detector ${detector} ${nonPrompt} ${combine} ${addOthers} " &

selection='lepSel3-onZ-njet3p-nbjet1p-Zpt0-leptonIso3'
nohup krenew -t -K 10 -- bash -c "python signalRegion.py ${small} --process ttZ --selection ${selection} --sample fwlite_ttZ_ll_LO_order2_15weights_ref --level ${level} --version ${version} --parameters ctZ 2.4 --luminosity 3000 --detector phase2_CMS --scale14TeV --addNonPrompt --nonInfoSignal" &
nohup krenew -t -K 10 -- bash -c "python signalRegion.py ${small} --process ttZ --selection ${selection} --sample fwlite_ttZ_ll_LO_order2_15weights_ref --level ${level} --version ${version} --parameters ctZI 2.4 --luminosity 3000 --detector phase2_CMS --scale14TeV --addNonPrompt --nonInfoSignal" &
nohup krenew -t -K 10 -- bash -c "python signalRegion.py ${small} --process ttZ --selection ${selection} --sample fwlite_ttZ_ll_LO_order2_15weights_ref --level ${level} --version ${version} --parameters cpQM 5 --luminosity 3000 --detector phase2_CMS --scale14TeV --addNonPrompt --nonInfoSignal" &
nohup krenew -t -K 10 -- bash -c "python signalRegion.py ${small} --process ttZ --selection ${selection} --sample fwlite_ttZ_ll_LO_order2_15weights_ref --level ${level} --version ${version} --parameters cpt 5 --luminosity 3000 --detector phase2_CMS --scale14TeV --addNonPrompt --nonInfoSignal" &
nohup krenew -t -K 10 -- bash -c "python signalRegion.py ${small} --process ttZ --selection ${selection} --sample fwlite_ttZ_ll_LO_order2_15weights_ref --level ${level} --version ${version} --luminosity 3000 --detector phase2_CMS --scale14TeV --addNonPrompt --nonInfoSignal" &


nohup krenew -t -K 10 -- bash -c "python signalRegion.py ${small} --process ttZ --selection ${selection} --sample fwlite_ttZ_ll_LO_order2_15weights_ref --level ${level} --version ${version} --luminosity 78 --detector phase2_CMS --addNonPrompt --combineTTX --addOthers --topEFTcolors" &

selection='lepSel3-onZ-njet3p-nRun2bjet1p-Zpt0-leptonIso3'
nohup krenew -t -K 10 -- bash -c "python signalRegion.py ${small} --process ttZ --selection ${selection} --sample fwlite_ttZ_ll_LO_order2_15weights_ref --level ${level} --version ${version} --luminosity 78 --detector CMS --addNonPrompt --combineTTX --addOthers --topEFTcolors" &
