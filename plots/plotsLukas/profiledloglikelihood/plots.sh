small=""
level="reco"
version="paper_240818"

selection='lepSel3-onZ-njet3p-nbjet1p-Zpt0-leptonIso3'
#nohup krenew -t -K 10 -- bash -c "python NLLPlot.py ${small} --process ttZ_3l --selection ${selection} --sample fwlite_ttZ_ll_LO_order2_15weights_ref --level ${level} --version ${version} --variables cpQM cpt --zRange 0 40 --luminosity 150 --binning 31 -8 28 31 -24 14 --fit BestFit --smooth --contours --cores 4" &
#nohup krenew -t -K 10 -- bash -c "python NLLPlot.py ${small} --process ttZ_3l --selection ${selection} --sample fwlite_ttZ_ll_LO_order2_15weights_ref --level ${level} --version paper_230818_zoom --variables cpQM cpt --zRange 0 40 --luminosity 150 --binning 31 -4 4 31 -4 4 --fit BestFit --smooth --contours --cores 4" &
#nohup krenew -t -K 10 -- bash -c "python NLLPlot.py ${small} --process ttZ_3l --selection ${selection} --sample fwlite_ttZ_ll_LO_order2_15weights_ref --level ${level} --version ${version} --variables ctZ ctZI --zRange 0 40 --luminosity 150 --binning 31 -1.5 1.5 31 -1.5 1.5 --fit BestFit --smooth --contours --cores 4" &

#CuB CuW
#nohup krenew -t -K 10 -- bash -c "python NLLPlot_CuB_CuW.py ${small} --process ttZ_3l --selection ${selection} --sample fwlite_ttZ_ll_LO_order2_15weights_ref --level ${level} --version ${version} --zRange 0 40 --luminosity 150 --binning 31 -2 2 31 -2 2 --fit BestFit --smooth --contours --cores 4" &


selection='lepSel1-gammapt40-njet3p-nbjet1p-relIso0to0.12-met40-leptonIso1'
#nohup krenew -t -K 10 -- bash -c "python NLLPlot.py ${small} --process ttgamma_1l --selection ${selection} --sample fwlite_ttgammaLarge_LO_order2_15weights_ref --level ${level} --version ${version} --variables ctZ ctZI --zRange 0 40 --luminosity 150 --binning 31 -0.5 0.5 31 -0.5 0.5 --smooth --contours --cores 4 --fit BestFit" &

#CuB CuW
#nohup krenew -t -K 10 -- bash -c "python NLLPlot_CuB_CuW.py ${small} --process ttgamma_1l --selection ${selection} --sample fwlite_ttgammaLarge_LO_order2_15weights_ref --level ${level} --version ${version} --zRange 0 40 --luminosity 150 --binning 31 -2 2 31 -2 2 --smooth --contours --cores 4 --fit BestFit" &


selection='lepSel2-gammapt40-njet2p-nbjet1p-relIso0to0.12-met40-leptonIso2'
#nohup krenew -t -K 10 -- bash -c "python NLLPlot.py ${small} --process ttgamma_2l --selection ${selection} --sample fwlite_ttgammaLarge_LO_order2_15weights_ref --level ${level} --version ${version} --variables ctZ ctZI --zRange 0 40 --luminosity 150 --binning 31 -0.5 0.5 31 -0.5 0.5 --smooth --contours --cores 4 --fit BestFit" &

#CuB CuW
#nohup krenew -t -K 10 -- bash -c "python NLLPlot_CuB_CuW.py ${small} --process ttgamma_2l --selection ${selection} --sample fwlite_ttgammaLarge_LO_order2_15weights_ref --level ${level} --version ${version} --zRange 0 40 --luminosity 150 --binning 31 -2 2 31 -2 2 --smooth --contours --cores 4 --fit BestFit" &


#combined
nohup krenew -t -K 10 -- bash -c "python NLLPlot_combined.py ${small} --level ${level} --version ${version} --variables ctZ ctZI --binning 31 -0.3 0.3 31 -0.3 0.3 --luminosity 150 --fit BestFit --zRange 0 40 --contours --smooth --cores 4" &
#python NLLPlot_combined.py ${small} --level ${level} --version ${version} --variables ctZ ctZI --binning 31 -0.3 0.3 31 -0.3 0.3 --luminosity 150 --fit BestFit --zRange 0 40 --contours --smooth --cores 4

#CuB CuW
nohup krenew -t -K 10 -- bash -c "python NLLPlot_combined_CuB_CuW.py ${small} --level ${level} --version ${version} --zRange 0 40 --luminosity 150 --binning 31 -2 2 31 -2 2 --smooth --contours --cores 4 --fit BestFit" &
#python NLLPlot_combined_CuB_CuW.py ${small} --level ${level} --version ${version} --zRange 0 40 --luminosity 150 --binning 31 -2 2 31 -2 2 --smooth --contours --cores 4 --fit BestFit
