version="approval"

selection='lepSel3-onZ-njet3p-nbjet1p-Zpt0-leptonIso3'
#nohup krenew -t -K 10 -- bash -c "python ptZ_uncertainties.py --selection ${selection} --version ${version} --parameters ctZ 2.0 ctZI 2.0 --shape --addUncertainties --scale14TeV --nonInfoSignal" &

echo "python ptZ_uncertainties.py --selection ${selection} --version ${version} --parameters ctZ 2.0 ctZI 2.0 --shape --addUncertainties --scale14TeV --nonInfoSignal" &

selection='lepSel3-onZ-njet3p-nbjet1p-Zpt200-leptonIso3'
#nohup krenew -t -K 10 -- bash -c "python cos_uncertainties.py --selection ${selection} --version ${version} --parameters ctZ 2.0 ctZI 2.0 --shape --addUncertainties --scale14TeV --nonInfoSignal" &

