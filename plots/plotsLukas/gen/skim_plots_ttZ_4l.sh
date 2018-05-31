#!/bin/bash

##############################################

#parameters
cpQM=$1
cpt=$2
ctW=$3
ctWI=$4
ctZ=$5
ctZI=$6
ctG=$7
ctGI=$8

# declare samples to analyze
#declare -a samples2=('fwlite_ttZ_ll_LO_order2_15weights')
declare -a samples2=('fwlite_ttZ_ll_LO_order2_15weights' 'fwlite_ttZ_ll_LO_order2_15weights_ref')
#declare -a samples2=('')
declare -a samples3=('fwlite_ttZ_ll_LO_order3_8weights')
#declare -a samples3=('')

# declare selection strings to analyze
declare -a selections=('lepSel4-onZ-njet2p-nbjet1p-Zpt0' 'lepSel4-onZ-njet2p-nbjet1p-Zpt0to100' 'lepSel4-onZ-njet2p-nbjet1p-Zpt100to200' 'lepSel4-onZ-njet2p-nbjet1p-Zpt200to300' 'lepSel4-onZ-njet2p-nbjet1p-Zpt300to400' 'lepSel4-onZ-njet2p-nbjet1p-Zpt400')

# declare sample size to analyze
declare -a samplesizes=('--small' '')
#declare -a samplesizes=('--small')
#declare -a samplesizes=('')

# define program to run by python
prog=skim_plots_ttZ_4l.py

#################################################

# 2nd order samples
order=2
for sample in "${samples2[@]}"
do

   if [ -z $sample ]; then
      continue
   fi

   for samplesize in "${samplesizes[@]}"
   do
      for selection in "${selections[@]}"
      do

         if [ -z $selection ]; then
            continue
         fi

         submitBatch.py --dpm "python ${prog} ${samplesize}             --sample ${sample} --order ${order} --selection ${selection} --parameters cpQM ${cpQM} cpt ${cpt} ctW ${ctW} ctWI ${ctWI} ctZ ${ctZ} ctZI ${ctZI} ctG ${ctG} ctGI ${ctGI}"
         submitBatch.py --dpm "python ${prog} ${samplesize} --scaleLumi --sample ${sample} --order ${order} --selection ${selection} --parameters cpQM ${cpQM} cpt ${cpt} ctW ${ctW} ctWI ${ctWI} ctZ ${ctZ} ctZI ${ctZI} ctG ${ctG} ctGI ${ctGI}"

      done
   done
done

# 3rd order samples
order=3
for sample in "${samples3[@]}"
do
   if [ -z $sample ]; then
      continue
   fi

   for samplesize in "${samplesizes[@]}"
   do
      for selection in "${selections[@]}"
      do

         if [ -z $selection ]; then
            continue
         fi

         submitBatch.py --dpm "python ${prog} ${samplesize}             --sample ${sample} --order ${order} --selection ${selection} --parameters cpQM ${cpQM} cpt ${cpt} ctW ${ctW} ctWI ${ctWI} ctZ ${ctZ} ctZI ${ctZI} ctG ${ctG} ctGI ${ctGI}"
         submitBatch.py --dpm "python ${prog} ${samplesize} --scaleLumi --sample ${sample} --order ${order} --selection ${selection} --parameters cpQM ${cpQM} cpt ${cpt} ctW ${ctW} ctWI ${ctWI} ctZ ${ctZ} ctZI ${ctZI} ctG ${ctG} ctGI ${ctGI}"

      done
   done
done