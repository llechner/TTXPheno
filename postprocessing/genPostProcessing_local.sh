#!/bin/sh

samplename=fwlite_ttZ_ll_LO_order2_15weights_ref
order=2
nJobs=20

for i in $(seq 0 $((${nJobs}-1)))
do
   python genPostProcessing.py --nJobs ${nJobs} --job ${i} --overwrite --logLevel DEBUG --sample ${samplename} --addReweights --interpolationOrder ${order} --delphes > log_${i}.log 2>&1 &
done


