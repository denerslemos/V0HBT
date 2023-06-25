# V0HBT

V0 HBT in pPb at 8 TeV
Setup CMSSW (work inside of CMSSW80X):
```
ssh -XY user@lxplus.cern.ch
export SCRAM_ARCH=slc7_amd64_gcc530
cmsrel CMSSW_8_0_24
cd CMSSW_8_0_24/src
cmsenv
```
Download: 
```
git clone https://github.com/denerslemos/V0HBT.git
```
Compile:
```
cd V0HBT
g++ hbtV0.C hbt.h `root-config --libs` `root-config --cflags` -o executable.exe
```
Submit jobs (first edit the script "submit_new.py" to have your desired configuration, see GDoc: https://docs.google.com/document/d/1E7WP24n_zG8nXQ49bEcL2uHEr0VMtwb29931xo5AEiA/edit):
```
python submit_new.py
```

After jobs are done, we merge all histograms using HADD (example for V0loose):
```
hadd Merged_V0HBT_AAAAA.root V0loose/*.root Pbp/V0loose/*.root 
```
where ```AAAAA``` corresponds to the systematics. We will need make sure to clean up the folder after merging histograms are big.

