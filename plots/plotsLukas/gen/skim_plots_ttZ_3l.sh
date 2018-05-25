#!/bin/bash

##############################################

#parameters
cpQM=10
cpt=10
ctW=10
ctWI=10
ctZ=10
ctZI=10
ctG=10
ctGI=10

# declare samples to analyze
#declare -a samples2=('fwlite_ttZ_ll_LO_order2_15weights')
declare -a samples2=('')
declare -a samples3=('fwlite_ttZ_ll_LO_order3_8weights')

# declare selection strings to analyze
declare -a selections=('lepSel3-onZ-njet3p-nbjet1p-Zpt0', 'lepSel3-onZ-njet3p-nbjet1p-Zpt100to200', 'lepSel3-onZ-njet3p-nbjet1p-Zpt200to300', 'lepSel3-onZ-njet3p-nbjet1p-Zpt300to400', 'lepSel3-onZ-njet3p-nbjet1p-Zpt400')

# declare sample size to analyze
declare -a samplesizes=('--small') #, '')

# declare scales to analyze
declare -a scales=('--scaleLumi', '')

# define program to run by python
prog=skim_plots_ttZ_3l.py

#################################################

# 2nd order samples
order=2
for sample in "${samples2[@]}"
do
   for samplesize in "${samplesizes[@]}"
   do
      for selection in "${selections[@]}"
      do
         for scale in "${scales[@]}"
         do
            python ${prog} ${samplesize} ${scale} --sample ${sample} --order ${order} --selection ${selection} --parameters cpQM ${cpQM} cpt ${cpt} ctW ${ctW} ctWI ${ctWI} ctZ ${ctZ} ctZI ${ctZI} ctG ${ctG} ctGI ${ctGI}
         done
      done
   done
done

# 3rd order samples
order=3
for sample in "${samples3[@]}"
do
   for samplesize in "${samplesizes[@]}"
   do
      for selection in "${selections[@]}"
      do
         for scale in "${scales[@]}"
         do
            python ${prog} ${samplesize} ${scale} --sample ${sample} --order ${order} --selection ${selection} --parameters cpQM ${cpQM} cpt ${cpt} ctW ${ctW} ctWI ${ctWI} ctZ ${ctZ} ctZI ${ctZI} ctG ${ctG} ctGI ${ctGI}
         done
      done
   done
done
