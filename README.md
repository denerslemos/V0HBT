# V0HBT

V0 HBT in pPb at 8 TeV
Setup CMSSW (work inside of CMSSW80X):
```
cmsrel CMSSW_8_0_24
cd CMSSW_8_0_24/src
cmsenv
```
Download: 
```
git clone git@github.com:denerslemos/V0HBT.git
```
Compile:
```
cd V0HBT
g++ hbtV0.C hbt.h `root-config --libs` `root-config --cflags` -o executable.exe
```
Submit jobs:
```
python submit_new.py
```

