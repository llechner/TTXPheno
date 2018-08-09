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
declare -a samples2=('fwlite_ttgamma_LO_order2_15weights_ref')
declare -a samples3=('')

# declare selection strings to analyze
#declare -a selections=('gammapt40-nlep1p-njet3p-nbjet1p' 'gammapt40to100-nlep1p-njet3p-nbjet1p' 'gammapt100to200-nlep1p-njet3p-nbjet1p' 'gammapt200to300-nlep1p-njet3p-nbjet1p' 'gammapt300-nlep1p-njet3p-nbjet1p')
#declare -a selections=('gammapt40-nlep1p-njet3p-nbjet1p' 'gammapt40-nlepttgamma2l1p' 'gammapt15-nlep1p-njet3p-nbjet1p')
#declare -a selections=('lepSel1-gammapt40-njet3p-nbjet1p-relIso0to0.1' 'lepSel1-gammapt15-njet2p-nbjet1-relIso0to0.1' 'lepSel1-gammapt15-njet2p-nbjet2p-relIso0to0.1' 'lepSel1-gammapt15-njet3p-nbjet1p-relIso0to0.1')
#declare -a selections=('lepSel1-gammapt40-njet3p-nbjet1p-relIso0to0.4' 'lepSel1-gammapt15-njet2p-nbjet1-relIso0to0.4' 'lepSel1-gammapt15-njet2p-nbjet2p-relIso0to0.4' 'lepSel1-gammapt15-njet3p-nbjet1p-relIso0to0.4')
#declare -a selections=('lepSel1-gammapt40-njet3p-nbjet1p-relIso0to0.4-met40' 'lepSel1-gammapt100-njet3p-nbjet1p-relIso0to0.4-met40' 'lepSel1-gammapt40-njet3p-nbjet1p-relIso0to0.12-met40' 'lepSel1-gammapt100-njet3p-nbjet1p-relIso0to0.12-met40')
#declare -a selections=('lepSel1-gammapt40-njet3p-nbjet1p-relIso0to0.4-met40')
#declare -a selections=('lepSel1-gammapt40-njet3p-nbjet1p-relIso0to0.12-met40' 'lepSel1-gammapt100-njet3p-nbjet1p-relIso0to0.12-met40' 'lepSel1-gammapt40-njet3p-nbjet1p-relIso0to0.4-met40' 'lepSel1-gammapt100-njet3p-nbjet1p-relIso0to0.4-met40')
declare -a selections=('lepSel1-gammapt200-njet3p-nbjet1p-relIso0to0.12-met40')

# declare sample size to analyze
#declare -a samplesizes=('--small' '')
#declare -a samplesizes=('--small')
declare -a samplesizes=('')

# declare reweighting
#declare -a reweightings=('' '--reweightPtXToSM')
#declare -a reweightings=('--reweightPtXToSM')
declare -a reweightings=('')

# declare scale
declare -a scales=('' '--scaleLumi')
#declare -a scales=('--scaleLumi')
#declare -a scales=('')

#declare -a levels=('gen')
#declare -a levels=('reco')
declare -a levels=('gen' 'reco')

#declare -a variables=("cpt" "cpQM")
declare -a variables=("cpt")

#declare -a flavors=('all' 'mu' 'e')
declare -a flavors=('all')

#declare -a binThresholds=("400" "100" "25" "0")
#declare -a binThresholds=("100" "0")
declare -a binThresholds=("100")
#declare -a binThresholds=("0")

#declare -a fisherInfo=("--addFisherInformation" "")
#declare -a fisherInfo=("--addFisherInformation")
declare -a fisherInfo=("")

#declare -a backgrounds=("--backgrounds" "")
declare -a backgrounds=("--backgrounds")

version='TTXPheno_08082018'
luminosity='150'
process='ttgamma_1l'

# define program to run by python
prog=skim_plots.py

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

            for level in "${levels[@]}"
            do

               for flavor in "${flavors[@]}"
               do

                  for variable in "${variables[@]}"
                  do

                     for binThreshold in "${binThresholds[@]}"
                     do

                        for addFisher in "${fisherInfo[@]}"
                        do

                           for background in "${backgrounds[@]}"
                           do

                               order=2
                               for sample in "${samples2[@]}"
                               do

                                  if [ -z $sample ]; then
                                     continue
                                  fi

#                                  echo "python ${prog} --processFile ${process} --luminosity ${luminosity} --version ${version} --level ${level} ${samplesize} ${reweight} ${scale} --sample ${sample} --order ${order} --selection ${selection} ${backgrounds} --parameters cpQM ${cpQM} cpt ${cpt} ctW ${ctW} ctWI ${ctWI} ctZ ${ctZ} ctZI ${ctZI} ctG ${ctG} ctGI ${ctGI} ${background} ${addFisher} --binThreshold ${binThreshold} --variables ${variable} --leptonFlavor ${flavor}"
                                  submitBatch.py --dpm "python ${prog} --processFile ${process} --luminosity ${luminosity} --version ${version} --level ${level} ${samplesize} ${reweight} ${scale} --sample ${sample} --order ${order} --selection ${selection} ${backgrounds} --parameters cpQM ${cpQM} cpt ${cpt} ctW ${ctW} ctWI ${ctWI} ctZ ${ctZ} ctZI ${ctZI} ctG ${ctG} ctGI ${ctGI} ${background} ${addFisher} --binThreshold ${binThreshold} --variables ${variable} --leptonFlavor ${flavor}"

                               done

                               order=3
                               for sample in "${samples3[@]}"
                               do

                                  if [ -z $sample ]; then
                                     continue
                                  fi

                                  submitBatch.py --dpm "python ${prog} --processFile ${process} --luminosity ${luminosity} --version ${version} --level ${level} ${samplesize} ${reweight} ${scale} --sample ${sample} --order ${order} --selection ${selection} ${backgrounds} --parameters cpQM ${cpQM} cpt ${cpt} ctW ${ctW} ctWI ${ctWI} ctZ ${ctZ} ctZI ${ctZI} ctG ${ctG} ctGI ${ctGI} ${background} ${addFisher} --binThreshold ${binThreshold} --variables ${variable} --leptonFlavor ${flavor}"

                               done
                            done
                        done
                     done
                  done
               done
            done
         done
      done
   done
done
