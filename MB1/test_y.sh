#!/bin/bash

echo "setup cmssw"
cd /afs/cern.ch/work/d/ddesouza/UIC/V0HBT/CMSSW_13_0_5/src/V0HBT/
cd /afs/cern.ch/work/d/ddesouza/UIC/V0HBT/CMSSW_13_0_5/src/
eval `scramv1 runtime -sh`
cd /afs/cern.ch/work/d/ddesouza/UIC/V0HBT/CMSSW_13_0_5/src/V0HBT//MB1
echo PWD: $PWD

./run.py -i $1