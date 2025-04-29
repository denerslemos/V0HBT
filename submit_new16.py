#!/usr/bin/env python

import os.path

V0home = "/afs/cern.ch/work/d/ddesouza/UIC/V0HBT/CMSSW_13_0_5/src/V0HBT/"
cmssource = "/afs/cern.ch/work/d/ddesouza/UIC/V0HBT/CMSSW_13_0_5/src/"

Nmix = 20  # number of events to mix (default: 20)
vzmix = 2.0  # vz matching in mixing (default: 2.0)
# mother removal (default: 3 - remove based on V0chi2; 0 - remove both; 1
# - remove randomly; 2 - remove based on mass; >3 - do not remove any
# mother)
method = 3
# systematic source (default: 1 - tight V0; 2 - loose V0; 9 - misid+; 10 -
# misid-; 14 - z vertex between |[3,15]|; 15 - zvertex [-3,3]; 16 - 1.5
# sigma peak ; 17 - 2.5 sigma peak; 18 - 3.5 sideband; 19 - 4.5 sideband;
# 23 - PU pileUpFilter_pPb8TeV_Gplus; 24 - PU pileUpFilter_pPb8TeV_vtx1;)
syst = 18
# true it submit jobs, false it just make the files (good for testing)
submit_jobs = True

folders = [
    'MB1',
    'MB2',
    'MB3',
    'MB4',
    'MB5',
    'MB6',
    'MB7',
    'MB8',
    'HM0120',
    'HM0150',
    'HM1',
    'HM2',
    'HM3',
    'HM4',
    'HM5',
    'HM6',
    'HM7',
    'Pbp/HM0120',
    'Pbp/HM0150',
    'Pbp/HM1',
    'Pbp/HM2',
    'Pbp/HM3',
    'Pbp/HM4',
    'Pbp/HM5',
    'Pbp/HM6',
    'Pbp/HM7',
    'Pbp/MB1',
    'Pbp/MB2',
    'Pbp/MB3',
    'Pbp/MB4',
    'Pbp/MB5',
    'Pbp/MB6',
    'Pbp/MB7',
    'Pbp/MB8',
    'Pbp/MB9',
    'Pbp/MB10',
    'Pbp/MB11',
    'Pbp/MB12',
    'Pbp/MB13',
    'Pbp/MB14',
    'Pbp/MB15',
    'Pbp/MB16',
    'Pbp/MB17',
    'Pbp/MB18',
    'Pbp/MB19',
    'Pbp/MB20']

os.system('chmod 777 executable.exe')
for i in folders:
    os.chdir(str(V0home) + '/' + i)
    os.system('mkdir -p V0std V0tight V0loose V0nhit3 V0pix0 V0dcaplus V0dcaplusminus V0misplus V0misminus V0miseeplus V0miseeminus V0rhoplus V0zvxt3to15 V0zvxt3 V0sigmapeakminus V0sigmapeakplus V0sigmasideplus V0sigmasideminus V0epos V0hijing V0ampt V0Gplus V0vtx1 V0noPU')
    os.system('rm V0std/* V0tight/* V0loose/* V0nhit3/* V0pix0/* V0dcaplus/* V0dcaplusminus/* V0misplus/* V0misminus/* V0miseeplus/* V0miseeminus/* V0rhoplus/* V0zvxt3to15/* V0zvxt3/* V0sigmapeakminus/* V0sigmapeakplus/* V0sigmasideplus/* V0sigmasideminus/* V0epos/* V0hijing/* V0ampt/* V0Gplus/* V0vtx1/* V0noPU/*')
    os.system('rm run*.py *.sub *.cc *.exe cond/*')
    os.system('cp ' + str(V0home) + '/executable.exe .')
    ff = open('test_y.sh', "w")
    shellfile = '''#!/bin/bash

echo "setup cmssw"
cd ''' + str(V0home) + '''
cd ''' + str(cmssource) + '''
eval `scramv1 runtime -sh`
cd ''' + str(V0home) + '''/''' + str(i) + '''
echo PWD: $PWD

./run.py -i $1'''
    ff.write(shellfile)
    ff.close()
    f = open('run.py', "w")
    command_lines = '''#!/usr/bin/env python3

from math import *
import os.path
import glob
import time
import optparse

Nmix = ''' + str(Nmix) + '''
vzmix = ''' + str(vzmix) + '''
meth = ''' + str(method) + '''
'''
    if i == 'MB1' or i == 'MB2' or i == 'MB3' or i == 'MB4' or i == 'MB5' or i == 'MB6' or i == 'MB7' or i == 'MB8' or i == 'Pbp/MB1' or i == 'Pbp/MB2' or i == 'Pbp/MB3' or i == 'Pbp/MB4' or i == 'Pbp/MB5' or i == 'Pbp/MB6' or i == 'Pbp/MB7' or i == 'Pbp/MB8' or i == 'Pbp/MB9' or i == 'Pbp/MB10' or i == 'Pbp/MB11' or i == 'Pbp/MB12' or i == 'Pbp/MB13' or i == 'Pbp/MB14' or i == 'Pbp/MB15' or i == 'Pbp/MB16' or i == 'Pbp/MB17' or i == 'Pbp/MB18' or i == 'Pbp/MB19' or i == 'Pbp/MB20':
        command_lines += '''mult = 0\n'''
    if i == 'HM0120' or i == 'Pbp/HM0120':
        command_lines += '''mult = 1\n'''
    if i == 'HM0150' or i == 'Pbp/HM0150':
        command_lines += '''mult = 2\n'''
    if i == 'HM1' or i == 'HM2' or i == 'HM3' or i == 'HM4' or i == 'HM5' or i == 'HM6' or i == 'Pbp/HM1' or i == 'Pbp/HM2' or i == 'Pbp/HM3' or i == 'Pbp/HM4' or i == 'Pbp/HM5' or i == 'Pbp/HM6':
        command_lines += '''mult = 3\n'''
    if i == 'HM7' or i == 'Pbp/HM7':
        command_lines += '''mult = 4\n'''
    if i == 'MB1' or i == 'MB2' or i == 'MB3' or i == 'MB4' or i == 'MB5' or i == 'MB6' or i == 'MB7' or i == 'MB8' or i == 'HM0120' or i == 'HM0150' or i == 'HM1' or i == 'HM2' or i == 'HM3' or i == 'HM4' or i == 'HM5' or i == 'HM6' or i == 'HM7':
        command_lines += '''beamside = 0\n'''
    else:
        command_lines += '''beamside = 1\n'''
#        command_lines += '''beamside = 2\n'''

    command_lines += '''syst = ''' + str(syst) + ''' #ignore 5 and 6
mcc = 0 #always 0 here -> only data

usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-i', '--inN', dest='inN', help='input number', default='0', type='int')
(opt, args) = parser.parse_args()
inn = opt.inN
start = time.time()\n'''
    if i == 'MB1':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/pPb/dataMB1pPb/pPb_MB1_'+str(inn)+'.root' \n'''
    if i == 'MB2':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/pPb/dataMB2pPb/pPb_MB2_'+str(inn)+'.root' \n'''
    if i == 'MB3':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/pPb/dataMB3pPb/pPb_MB3_'+str(inn)+'.root' \n'''
    if i == 'MB4':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/pPb/dataMB4pPb/pPb_MB4_'+str(inn)+'.root' \n'''
    if i == 'MB5':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/pPb/dataMB5pPb/pPb_MB5_'+str(inn)+'.root' \n'''
    if i == 'MB6':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/pPb/dataMB6pPb/pPb_MB6_'+str(inn)+'.root' \n'''
    if i == 'MB7':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/pPb/dataMB7pPb/pPb_MB7_'+str(inn)+'.root' \n'''
    if i == 'MB8':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/pPb/dataMB8pPb/pPb_MB8_'+str(inn)+'.root' \n'''
    if i == 'HM0120':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/Pbp/dataHM0120pPb/pPb_HM0120_'+str(inn)+'.root' \n'''
    if i == 'HM0150':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/Pbp/dataHM0150pPb/pPb_HM0150_'+str(inn)+'.root' \n'''
    if i == 'HM1':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/Pbp/dataHM1pPb/pPb_HM1_'+str(inn)+'.root' \n'''
    if i == 'HM2':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/Pbp/dataHM2pPb/pPb_HM2_'+str(inn)+'.root' \n'''
    if i == 'HM3':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/Pbp/dataHM3pPb/pPb_HM3_'+str(inn)+'.root' \n'''
    if i == 'HM4':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/Pbp/dataHM4pPb/pPb_HM4_'+str(inn)+'.root' \n'''
    if i == 'HM5':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/Pbp/dataHM5pPb/pPb_HM5_'+str(inn)+'.root' \n'''
    if i == 'HM6':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/Pbp/dataHM6pPb/pPb_HM6_'+str(inn)+'.root' \n'''
    if i == 'HM7':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/Pbp/dataHM7pPb/pPb_HM7_'+str(inn)+'.root' \n'''
    if i == 'Pbp/MB1':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/Pbp/dataMB1Pbp/Pbp_MB1_'+str(inn)+'.root' \n'''
    if i == 'Pbp/MB2':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/Pbp/dataMB2Pbp/Pbp_MB2_'+str(inn)+'.root' \n'''
    if i == 'Pbp/MB3':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/Pbp/dataMB3Pbp/Pbp_MB3_'+str(inn)+'.root' \n'''
    if i == 'Pbp/MB4':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/Pbp/dataMB4Pbp/Pbp_MB4_'+str(inn)+'.root' \n'''
    if i == 'Pbp/MB5':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/Pbp/dataMB5Pbp/Pbp_MB5_'+str(inn)+'.root' \n'''
    if i == 'Pbp/MB6':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/Pbp/dataMB6Pbp/Pbp_MB6_'+str(inn)+'.root' \n'''
    if i == 'Pbp/MB7':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/Pbp/dataMB7Pbp/Pbp_MB7_'+str(inn)+'.root' \n'''
    if i == 'Pbp/MB8':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/Pbp/dataMB8Pbp/Pbp_MB8_'+str(inn)+'.root' \n'''
    if i == 'Pbp/MB9':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/Pbp/dataMB9Pbp/Pbp_MB9_'+str(inn)+'.root' \n'''
    if i == 'Pbp/MB10':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/Pbp/dataMB10Pbp/Pbp_MB10_'+str(inn)+'.root' \n'''
    if i == 'Pbp/MB11':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/Pbp/dataMB11Pbp/Pbp_MB11_'+str(inn)+'.root' \n'''
    if i == 'Pbp/MB12':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/Pbp/dataMB12Pbp/Pbp_MB12_'+str(inn)+'.root' \n'''
    if i == 'Pbp/MB13':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/Pbp/dataMB13Pbp/Pbp_MB13_'+str(inn)+'.root' \n'''
    if i == 'Pbp/MB14':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/Pbp/dataMB14Pbp/Pbp_MB14_'+str(inn)+'.root' \n'''
    if i == 'Pbp/MB15':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/Pbp/dataMB15Pbp/Pbp_MB15_'+str(inn)+'.root' \n'''
    if i == 'Pbp/MB16':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/Pbp/dataMB16Pbp/Pbp_MB16_'+str(inn)+'.root' \n'''
    if i == 'Pbp/MB17':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/Pbp/dataMB17Pbp/Pbp_MB17_'+str(inn)+'.root' \n'''
    if i == 'Pbp/MB18':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/Pbp/dataMB18Pbp/Pbp_MB18_'+str(inn)+'.root' \n'''
    if i == 'Pbp/MB19':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/Pbp/dataMB19Pbp/Pbp_MB19_'+str(inn)+'.root' \n'''
    if i == 'Pbp/MB20':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/Pbp/dataMB20Pbp/Pbp_MB20_'+str(inn)+'.root' \n'''
    if i == 'Pbp/HM0120':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/Pbp/dataHM0120Pbp/Pbp_HM0120_'+str(inn)+'.root' \n'''
    if i == 'Pbp/HM0150':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/Pbp/dataHM0150Pbp/Pbp_HM0150_'+str(inn)+'.root' \n'''
    if i == 'Pbp/HM1':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/Pbp/dataHM1Pbp/Pbp_HM1_'+str(inn)+'.root' \n'''
    if i == 'Pbp/HM2':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/Pbp/dataHM2Pbp/Pbp_HM2_'+str(inn)+'.root' \n'''
    if i == 'Pbp/HM3':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/Pbp/dataHM3Pbp/Pbp_HM3_'+str(inn)+'.root' \n'''
    if i == 'Pbp/HM4':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/Pbp/dataHM4Pbp/Pbp_HM4_'+str(inn)+'.root' \n'''
    if i == 'Pbp/HM5':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/Pbp/dataHM5Pbp/Pbp_HM5_'+str(inn)+'.root' \n'''
    if i == 'Pbp/HM6':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/Pbp/dataHM6Pbp/Pbp_HM6_'+str(inn)+'.root' \n'''
    if i == 'Pbp/HM7':
        command_lines += '''file='/eos/cms/store/group/phys_heavyions/ddesouza/DATA/PUstd/Pbp/dataHM7Pbp/Pbp_HM7_'+str(inn)+'.root' \n'''

    command_lines += '''os.system('./executable.exe '+file+' '+str(Nmix)+' '+str(vzmix)+' '+str(meth)+' '+str(mult)+' '+str(inn)+' '+str(syst)+' '+str(mcc)+' '+str(beamside))

end =  time.time()
sec = end-start
a = sec//60//60//24
b = (sec//60//60)%24
c = (sec//60)%60
d = sec%60
print ("Execution Time:   "+str(int(a))+" days   "+str(int(b))+" hours   "+str(int(c))+" minutes   "+str(d)+" seconds")'''
    f.write(command_lines)
    f.close()
    os.system('chmod 777 run.py')
    if (submit_jobs):
        os.system('python3 test_y.py')
