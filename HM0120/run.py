#!/usr/bin/env python


from math import *
import os.path
import glob
import time
import optparse


Nmix = 20
vzmix = 2.0
meth = 3
mult = 1
beamside = 0
syst = 3 #ignore 5 and 6
mcc = 0 #always 0 here -> only data

usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-i', '--inN', dest='inN', help='input number', default='0', type='int')
(opt, args) = parser.parse_args()
inn = opt.inN
start = time.time()
file='/eos/user/d/ddesouza/DATA/PUstd/pPb/dataHM0120pPb/pPb_HM0120_'+str(inn)+'.root' 
os.system('./executable.exe '+file+' '+str(Nmix)+' '+str(vzmix)+' '+str(meth)+' '+str(mult)+' '+str(inn)+' '+str(syst)+' '+str(mcc)+' '+str(beamside))

end =  time.time()
sec = end-start
a = sec//60//60//24
b = (sec//60//60)%24
c = (sec//60)%60
d = sec%60
print "Execution Time:   "+str(int(a))+" days   "+str(int(b))+" hours   "+str(int(c))+" minutes   "+str(d)+" seconds"