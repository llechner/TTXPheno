#small="--small"
small=""
level="reco"
version="approval"

scale="--scale14TeV"
detector="phase2_CMS"

nonPrompt="--addNonPrompt"
#nonPrompt=""

#combine="--combineTTX"
combine=""

#addOthers="--addOthers"
addOthers=""

selection='lepSel3-onZ-njet3p-nbjet1p-Zpt0-leptonIso3'
nohup krenew -t -K 10 -- bash -c "python signalRegion.py ${small} --process ttZ --selection ${selection} --sample fwlite_ttZ_ll_LO_order2_15weights_ref --level ${level} --version ${version} --parameters ctZ 2 --luminosity 3000 --detector phase2_CMS --scale14TeV --addNonPrompt --nonInfoSignal" &
