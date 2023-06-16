#!/bin/bash

echo "setup cmssw"
cd /afs/cern.ch/work/d/ddesouza/ana
cd /afs/cern.ch/work/d/ddesouza/CMS_V0/Full/CMSSW_8_0_24/src
eval `scramv1 runtime -sh`
cd /afs/cern.ch/work/d/ddesouza/ana/Pbp/MB14
echo PWD: $PWD

./run.py -i $1