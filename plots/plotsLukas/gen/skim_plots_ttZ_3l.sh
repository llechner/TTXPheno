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
#declare -a samples2=('fwlite_ttZ_ll_LO_order2_15weights' 'fwlite_ttZ_ll_LO_order2_15weights_ref')
declare -a samples2=('fwlite_ttZ_ll_LO_order2_15weights_ref')
#declare -a samples2=('fwlite_ttZ_ll_LO_order2_15weights')
#declare -a samples2=('')

#declare -a samples3=('fwlite_ttZ_ll_LO_order3_8weights')
declare -a samples3=('')

# declare selection strings to analyze
declare -a selections=('lepSel3-onZ-njet3p-nbjet1p-Zpt0' 'lepSel3-onZ-njet3p-nbjet1p-Zpt0to100' 'lepSel3-onZ-njet3p-nbjet1p-Zpt100to200' 'lepSel3-onZ-njet3p-nbjet1p-Zpt200to300' 'lepSel3-onZ-njet3p-nbjet1p-Zpt300to400' 'lepSel3-onZ-njet3p-nbjet1p-Zpt400')

# declare sample size to analyze
#declare -a samplesizes=('--small' '')
#declare -a samplesizes=('--small')
declare -a samplesizes=('')

# declare reweighting
declare -a reweightings=('' '--reweightPtZToSM')
#declare -a reweightings=('--reweightPtZToSM')
#declare -a reweightings=('')

# declare scale
declare -a scales=('' '--scaleLumi')
#declare -a scales=('--scaleLumi')
#declare -a scales=('')

# define program to run by python
prog=skim_plots_ttZ_3l.py

#################################################

for samplesize in "${samplesizes[@]}"
do
   for selection in "${selections[@]}"
   do

      if [ -z $selection ]; then
         continue
      fi

      for scale in "${scales[@]}"
      do

         for reweight in "${reweightings[@]}"
         do

            order=2
            for sample in "${samples2[@]}"
            do

               if [ -z $sample ]; then
                  continue
               fi

               submitBatch.py --dpm "python ${prog} ${samplesize} ${reweight} ${scale} --sample ${sample} --order ${order} --selection ${selection} --parameters cpQM ${cpQM} cpt ${cpt} ctW ${ctW} ctWI ${ctWI} ctZ ${ctZ} ctZI ${ctZI} ctG ${ctG} ctGI ${ctGI}"

            done

            order=3
            for sample in "${samples3[@]}"
            do

               if [ -z $sample ]; then
                  continue
               fi

               submitBatch.py --dpm "python ${prog} ${samplesize} ${reweight} ${scale} --sample ${sample} --order ${order} --selection ${selection} --parameters cpQM ${cpQM} cpt ${cpt} ctW ${ctW} ctWI ${ctWI} ctZ ${ctZ} ctZI ${ctZI} ctG ${ctG} ctGI ${ctGI}"

            done
         done
      done
   done
done
