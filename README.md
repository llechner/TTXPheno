# TTXPheno
## Get code for TTXPheno studies

```
cmsrel CMSSW_9_4_3
cd CMSSW_9_4_3/src
cmsenv
git cms-init
git clone https://github.com/TTXPheno/TTXPheno
# looping & plotting - needed for the examples
git clone https://github.com/TTXPheno/RootTools
scram b -j40
```

## Prerequisites
Add your user-specific locations to [user.py](https://github.com/TTXPheno/TTXPheno/blob/master/Tools/python/user.py).
You can do that also later and thus fix potential 'ImportError' from `user.py`  when they occur. 

## How to make a flat ntuple from a GEN file
Define the sample in e.g. [fwlite_benchmarks.py](https://github.com/TTXPheno/TTXPheno/blob/master/samples/python/fwlite_benchmarks.py) and do 
```
voms-proxy-init -voms cms -out ~/private/.proxy
cd TTXPheno/postprocessing
python genPostProcessing.py --overwrite --sample fwlite_ttZ_ll_LO_highStat_scan --addReweights --small
```
Remove `--small` to process the full sample. Add `--nJobs=10` and `--job=0` etc. options to run only on the first 10% of events. In Vienna you could also do 
```
export X509_USER_PROXY=~/private/.proxy
submitBatch.py --dpm genPostProcessing.sh # append `#SPLIT200` to a line
``` 
For other places we can add submission scripts.

## Examples
