#!/bin/sh

#####################################################################################

# create bin and cut plots with plotting variables given in process_variables_ROC.py
# script loops over all indices in the range of [startIndex, endIndex]
# indices which are not specified are skipped in the python script

#small="--small"
small=""


#declare -a variables=("cpt" "cpQM" "ctZ" "ctW" "cpt cpQM" "ctZ ctZI ctW ctWI" "cpt cpQM ctZ ctZI" "cpt cpQM ctW ctWI" "ctZ ctW" "ctZ ctZI" "ctW ctWI" "cpt cpQM ctZ ctW")
#declare -a variables=("cpt cpQM ctZ ctZI" "cpt cpQM" "ctZ ctZI" "cpt cpQM ctZ ctW" "cpt cpQM ctW ctWI")
#declare -a variables=("ctZ ctZI" "ctZ ctW")
#declare -a variables=("cpt" "cpQM" "cpt cpQM" "cpt cpQM ctZ ctZI ctW ctWI")
#declare -a variables=("cpt" "cpQM" "cpt cpQM")
declare -a variables=("cpt cpQM")

#declare -a levels=("reco" "gen")
#declare -a levels=("reco" "genLep" "gen")
#declare -a levels=("gen")
declare -a levels=("reco")

#declare -a fpsScales=("--fpsScaling" "")
declare -a fpsScales=("")
#declare -a fpsScales=("--fpsScaling")

version="TTXPheno_08082018"

#declare -a binThresholds=("400" "100" "25" "0")
declare -a binThresholds=("100")

#declare -a parameters=("" "--parameters ctZ 3 ctZI 3")
declare -a parameters=("")


#####################################################################################

for variable in "${variables[@]}"
do

    for level in "${levels[@]}"
    do

        for fpsScale in "${fpsScales[@]}"
        do

            for binThreshold in "${binThresholds[@]}"
            do

                for parameter in "${parameters[@]}"
                do

#                    echo "python fisher_information_ev.py ${small} --level ${level} --version ${version} --sample fwlite_ttZ_ll_LO_order2_15weights_ref  --process ttZ     --order 2 --selection lepSel3-onZ-njet3p-nbjet1p-Zpt0 --variables ${variable} ${parameter} ${fpsScale} --binThreshold ${binThreshold}"
                    submitBatch.py --dpm "python fisher_information_ev.py ${small} --level ${level} --version ${version} --sample fwlite_ttZ_ll_LO_order2_15weights_ref  --process ttZ     --order 2 --selection lepSel3-onZ-njet3p-nbjet1p-Zpt0 --variables ${variable} ${parameter} ${fpsScale} --binThreshold ${binThreshold}"
                    submitBatch.py --dpm "python fisher_information_ev.py ${small} --level ${level} --version ${version} --sample fwlite_ttZ_ll_LO_order2_15weights_ref  --process ttZ     --order 2 --selection lepSel3-onZ-njet3p-nbjet1p-Zpt200 --variables ${variable} ${parameter} ${fpsScale} --binThreshold ${binThreshold}"
#                    submitBatch.py --dpm "python fisher_information_ev.py ${small} --level ${level} --version ${version} --sample fwlite_ttZ_ll_LO_order2_15weights_ref  --process ttZ     --order 2 --selection lepSel4-onZ-njet2p-nbjet1p-Zpt0 --variables ${variable} ${parameter} ${fpsScale} --binThreshold ${binThreshold}"
#                    submitBatch.py --dpm "python fisher_information_ev.py ${small} --level ${level} --version ${version} --sample fwlite_ttW_LO_order2_15weights_ref     --process ttW     --order 2 --selection nlep2p-njet2p-nbjet1p-Wpt0 --variables ${variable} ${parameter} ${fpsScale} --binThreshold ${binThreshold}"
#                    submitBatch.py --dpm "python fisher_information_ev.py ${small} --level ${level} --version ${version} --sample fwlite_ttgamma_LO_order2_15weights_ref --process ttgamma --order 2 --selection lepSel1-gammapt40-njet3p-nbjet1p-relIso0to0.4-met40 --variables ${variable} ${parameter} ${fpsScale} --binThreshold ${binThreshold}"
#                    submitBatch.py --dpm "python fisher_information_ev.py ${small} --level ${level} --version ${version} --sample fwlite_ttgamma_LO_order2_15weights_ref --process ttgamma --order 2 --selection lepSel2-gammapt40-njet2p-nbjet1p-relIso0to0.4 --variables ${variable} ${parameter} ${fpsScale} --binThreshold ${binThreshold}"

                done
            done
        done
    done
done
