# V0HBT

V0 HBT in pPb at 8 TeV
Setup CMSSW (work inside of CMSSW80X):
```
cmsrel CMSSW_13_0_5
cd CMSSW_13_0_5/src/
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
Go to the ```submit_newX.py``` and change the paths in lines 5 and 6. Then submit jobs:
```
python3 submit_newX.py
```
where "X" is changes from 0 to 19 for each systematic source. See the GDOCS.

After jobs are done, we merge all histograms using HADD (example for V0loose):
```
hadd Merged_V0HBT_AAAAA.root V0loose/*.root Pbp/V0loose/*.root 
```
where ```AAAAA``` corresponds to the systematics. We will need make sure to clean up the folder after merging histograms are big.


