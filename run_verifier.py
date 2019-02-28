import sys
import os
import tempfile
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--samplerType', type=int, help="Sampler Type;\n UniGen=1, Solver=2, QuickSampler=3,STS=4, OtherSampler=5;\n" + \
  "(Please make appropriate changes to code for OtherSampler);\n" + "default = 2", default = 2, dest= 'samplerType')
tempDir = tempfile.gettempdir()
processrankStr = os.environ.get("OMPI_COMM_WORLD_RANK")
if (processrankStr == None):
    processrankStr = 0
processrank = int(processrankStr)
#print(processrank)
#print(processrankStr)
all_file_name = str(tempDir)+'/all_files_'+str(processrank)+'.log'

cmd = 'find . -name \*.cnf >' + all_file_name

os.system(cmd)
f = open(all_file_name,'r')
lines = f.readlines()
f.close()
#samplerType = 4-(processrank/len(lines))
args = parser.parse_args()
samplerType = args.samplerType
if (samplerType < 1):
    exit(-1)
if (samplerType > 5):
    exit(-1)
filepos = lines[processrank%len(lines)].strip()
fileSuffix = filepos.split('/')[-1][:-4]
cmd ='mkdir outDir'
os.system(cmd)
#cmd = 'ulimit -v 4000000; python Verifier.py --sampler '+str(samplerType)+' --inverted 1 --reverse 0 --seed 0 --exp 1000 '+filepos+'  outDir/sampler_'+str(samplerType)+'_'+fileSuffix+'.out'
#print(cmd)
cmd = 'python Verifier.py --sampler '+str(samplerType)+' --reverse 0 --seed 0 --exp 1000 '+filepos+'  outDir/sampler_'+str(samplerType)+'_'+fileSuffix+'.out'
os.system(cmd)
