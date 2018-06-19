#!/bin/sh

#####################################################################################

# create bin and cut plots with plotting variables given in process_variables_ROC.py
# script loops over all indices in the range of [startIndex, endIndex]
# indices which are not specified are skipped in the python script

startIndex=1
endIndex=10

small='--small'
#small=''

declare -a variables=('cpt cpQM' 'ctZ ctZI' 'cpt cpQM ctZ ctZI ctW ctWI' )

#################################################

#####################################################################################

for variable in "${variables[@]}"
do

    submitBatch.py --dpm 'python fisher_information_full.py $small --sample fwlite_ttZ_ll_LO_order2_15weights_ref  --process ttZ     --order 2 --selection lepSel3-onZ-njet3p-nbjet1p-Zpt0 --variables cpt cpQM ctZ ctZI ctW ctWI'
    submitBatch.py --dpm 'python fisher_information_full.py $small --sample fwlite_ttW_LO_order2_15weights_ref     --process ttW     --order 2 --selection nlep2p-njet2p-nbjet1p-Wpt0      --variables cpt cpQM ctZ ctZI ctW ctWI'
    submitBatch.py --dpm 'python fisher_information_full.py $small --sample fwlite_ttgamma_LO_order2_15weights_ref --process ttgamma --order 2 --selection gammapt40-nlep1p-njet3p-nbjet1p --variables cpt cpQM ctZ ctZI ctW ctWI'

    submitBatch.py --dpm 'python fisher_information_fps.py $small --sample fwlite_ttZ_ll_LO_order2_15weights_ref  --process ttZ     --order 2 --selection lepSel3-onZ-njet3p-nbjet1p-Zpt0 --variables cpt cpQM ctZ ctZI ctW ctWI'
    submitBatch.py --dpm 'python fisher_information_fps.py $small --sample fwlite_ttW_LO_order2_15weights_ref     --process ttW     --order 2 --selection nlep2p-njet2p-nbjet1p-Wpt0      --variables cpt cpQM ctZ ctZI ctW ctWI'
    submitBatch.py --dpm 'python fisher_information_fps.py $small --sample fwlite_ttgamma_LO_order2_15weights_ref --process ttgamma --order 2 --selection gammapt40-nlep1p-njet3p-nbjet1p --variables cpt cpQM ctZ ctZI ctW ctWI'

done

for i in $(seq $startIndex $endIndex)
do

    submitBatch.py --dpm 'python fisher_information_ROC_binPlots.py $small --sample fwlite_ttZ_ll_LO_order2_15weights_ref  --process ttZ     --order 2 --selection lepSel3-onZ-njet3p-nbjet1p-Zpt0 --variables cpt cpQM ctZ ctZI ctW ctWI --selectPlots $i'
    submitBatch.py --dpm 'python fisher_information_ROC_cutPlots.py $small --sample fwlite_ttZ_ll_LO_order2_15weights_ref  --process ttZ     --order 2 --selection lepSel3-onZ-njet3p-nbjet1p-Zpt0 --variables cpt cpQM ctZ ctZI ctW ctWI --selectPlots $i'

    submitBatch.py --dpm 'python fisher_information_ROC_binPlots.py $small --sample fwlite_ttW_LO_order2_15weights_ref     --process ttW     --order 2 --selection nlep2p-njet2p-nbjet1p-Wpt0      --variables cpt cpQM ctZ ctZI ctW ctWI --selectPlots $i'
    submitBatch.py --dpm 'python fisher_information_ROC_cutPlots.py $small --sample fwlite_ttW_LO_order2_15weights_ref     --process ttW     --order 2 --selection nlep2p-njet2p-nbjet1p-Wpt0      --variables cpt cpQM ctZ ctZI ctW ctWI --selectPlots $i'

    submitBatch.py --dpm 'python fisher_information_ROC_binPlots.py $small --sample fwlite_ttgamma_LO_order2_15weights_ref --process ttgamma --order 2 --selection gammapt40-nlep1p-njet3p-nbjet1p --variables cpt cpQM ctZ ctZI ctW ctWI --selectPlots $i'
    submitBatch.py --dpm 'python fisher_information_ROC_cutPlots.py $small --sample fwlite_ttgamma_LO_order2_15weights_ref --process ttgamma --order 2 --selection gammapt40-nlep1p-njet3p-nbjet1p --variables cpt cpQM ctZ ctZI ctW ctWI --selectPlots $i'

done
