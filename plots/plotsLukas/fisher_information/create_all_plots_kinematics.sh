#!/bin/sh

#####################################################################################

# create bin and cut plots with plotting variables given in process_variables_ROC.py
# script loops over all indices in the range of [startIndex, endIndex]
# indices which are not specified are skipped in the python script

startIndex=1
endIndex=10

#small="--small"
small=""

declare -a variables=("ctZ" "ctW" "ctW ctWI" "cpt cpQM" "ctZ ctZI")
#declare -a variables=("ctWI ctZI" "ctW ctZ")

declare -a parameters=("ctZ 3 ctZI 3 ctW 3 ctWI 3")

declare -a levels=("reco" "genLep" "gen")

version="v7"

#################################################

#####################################################################################

for variable in "${variables[@]}"
do

    for level in "${levels[@]}"
    do

        for i in $(seq ${startIndex} ${endIndex})
        do

            submitBatch.py --dpm "python fisher_information_kinematics.py ${small} --level ${level} --version ${version} --sample fwlite_ttZ_ll_LO_order2_15weights_ref  --process ttZ     --order 2 --selection lepSel3-onZ-njet3p-nbjet1p-Zpt0 --variables ${variable} --selectPlots ${i}"
            submitBatch.py --dpm "python fisher_information_kinematics.py ${small} --level ${level} --version ${version} --sample fwlite_ttZ_ll_LO_order2_15weights_ref  --process ttZ     --order 2 --selection lepSel3-onZ-njet3p-nbjet1p-Zpt0 --variables ${variable} --selectPlots ${i}"

#            submitBatch.py --dpm "python fisher_information_kinematics.py ${small} --level ${level} --version ${version} --sample fwlite_ttW_LO_order2_15weights_ref     --process ttW     --order 2 --selection nlep2p-njet2p-nbjet1p-Wpt0      --variables ${variable} --selectPlots ${i}"
#            submitBatch.py --dpm "python fisher_information_kinematics.py ${small} --level ${level} --version ${version} --sample fwlite_ttW_LO_order2_15weights_ref     --process ttW     --order 2 --selection nlep2p-njet2p-nbjet1p-Wpt0      --variables ${variable} --selectPlots ${i}"

#            submitBatch.py --dpm "python fisher_information_kinematics.py ${small} --level ${level} --version ${version} --sample fwlite_ttgamma_LO_order2_15weights_ref --process ttgamma --order 2 --selection gammapt40-nlep1p-njet3p-nbjet1p --variables ${variable} --selectPlots ${i}"
#            submitBatch.py --dpm "python fisher_information_kinematics.py ${small} --level ${level} --version ${version} --sample fwlite_ttgamma_LO_order2_15weights_ref --process ttgamma --order 2 --selection gammapt40-nlep1p-njet3p-nbjet1p --variables ${variable} --selectPlots ${i}"

            for parameter in "${parameters[@]}"
            do

                submitBatch.py --dpm "python fisher_information_kinematics.py ${small} --level ${level} --version ${version} --sample fwlite_ttZ_ll_LO_order2_15weights_ref  --process ttZ     --order 2 --selection lepSel3-onZ-njet3p-nbjet1p-Zpt0 --variables ${variable} --selectPlots ${i} --parameters ${parameter}"
                submitBatch.py --dpm "python fisher_information_kinematics.py ${small} --level ${level} --version ${version} --sample fwlite_ttZ_ll_LO_order2_15weights_ref  --process ttZ     --order 2 --selection lepSel3-onZ-njet3p-nbjet1p-Zpt0 --variables ${variable} --selectPlots ${i} --parameters ${parameter}"

#                submitBatch.py --dpm "python fisher_information_kinematics.py ${small} --level ${level} --version ${version} --sample fwlite_ttW_LO_order2_15weights_ref     --process ttW     --order 2 --selection nlep2p-njet2p-nbjet1p-Wpt0      --variables ${variable} --selectPlots ${i} --parameters ${parameter}"
#                submitBatch.py --dpm "python fisher_information_kinematics.py ${small} --level ${level} --version ${version} --sample fwlite_ttW_LO_order2_15weights_ref     --process ttW     --order 2 --selection nlep2p-njet2p-nbjet1p-Wpt0      --variables ${variable} --selectPlots ${i} --parameters ${parameter}"

#                submitBatch.py --dpm "python fisher_information_kinematics.py ${small} --level ${level} --version ${version} --sample fwlite_ttgamma_LO_order2_15weights_ref --process ttgamma --order 2 --selection gammapt40-nlep1p-njet3p-nbjet1p --variables ${variable} --selectPlots ${i} --parameters ${parameter}"
#                submitBatch.py --dpm "python fisher_information_kinematics.py ${small} --level ${level} --version ${version} --sample fwlite_ttgamma_LO_order2_15weights_ref --process ttgamma --order 2 --selection gammapt40-nlep1p-njet3p-nbjet1p --variables ${variable} --selectPlots ${i} --parameters ${parameter}"

            done
        done
    done
done

