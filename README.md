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
Add your user-specific locations to `[[https://github.com/TTXPheno/TTXPheno/Tools/python/user.py][https://github.com/TTXPheno/TTXPheno/Tools/python/user.py]]`
Depending on what you're planning, you need not do that for all of locations and 
fix potential 'ImportError' from `TXPheno/Tools/python/user.py`  when they occur. 

## How to make a flat ntuple from a GEN file
