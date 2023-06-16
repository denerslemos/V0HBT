#!/bin/bash

echo "setup cmssw"
cd /afs/cern.ch/work/d/ddesouza/
#export SCRAM_ARCH=slc7_amd64_gcc700
cd /afs/cern.ch/work/d/ddesouza/CMS_V0/Full/CMSSW_8_0_24/src
#cmsenv
eval `scramv1 runtime -sh`
cd /afs/cern.ch/work/d/ddesouza/ana/HM5/
#g++ hbtV0.C hbt.h `root-config --libs` `root-config --cflags` -o executable.exe
echo PWD: $PWD

./run.py -i $1 

