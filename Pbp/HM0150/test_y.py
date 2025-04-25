#!/usr/bin/env python3

import subprocess
name = "submit_test.sub"
f = open(name, "w")

command_lines = '''universe   = vanilla
getenv     = True
executable = test_y.sh
+JobFlavour           = "tomorrow"
requirements = ((OpSysAndVer =?= "AlmaLinux9") && (CERNEnvironment =?= "qa"))
RequestCpus = 2
'''

i = 'HM' #sample folder
a = 20 #number of files in this folder
for j in range(0, a):
       temp = '''
log        = cond/test.'''+str(j+1)+'''.'''+str(i)+'''.log
output     = cond/test.'''+str(j+1)+'''.'''+str(i)+'''.out
error      = cond/test.'''+str(j+1)+'''.'''+str(i)+'''.err
arguments = '''+str(j+1)+'''
queue
   '''
       command_lines += temp

f.write(command_lines)
f.close()
subprocess.call(["condor_submit", name])
