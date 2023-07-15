import os
import sys
import numpy as np
from gsd import hoomd as gsd
import argparse
import os.path
import glob
import hoomd
import hoomd.md

parser=argparse.ArgumentParser()
parser.add_argument("--fileprefix",type=str)
parser.add_argument("--vfr",type=float)

args = parser.parse_args()
locals().update(vars(args))

#outputfile=fileprefix+".allruns.gsd"

hoomd.context.initialize("")
#traj = gsd.open(name=outputfile, mode='wb')

print("/"*80)
print("Fileprefix: ",os.path.basename(fileprefix))

k=0

if(vfr<=0.35):
    num_files=3
    for filename in sorted(glob.glob(fileprefix+"_N*.gsd"),key=lambda x:int(os.path.basename(x).split("_")[15].split('.')[0].replace('N',''))):
        if(k>num_files-1):
            break
        trajectory = gsd.open(filename,'rb') # read gsd file
        key=int(os.path.basename(filename).split("_")[15].split('.')[0].replace('N',''))
        print("Run no: ",key)
        print("No of frames: ",len(trajectory))
        """
        for i in range(len(trajectory)):
            if(k!=0 and i==0):
                continue
            system= trajectory[int(i)]
            traj.append(system) 
        k+=1
        """
else:
    #num_files=3
    for filename in sorted(glob.glob(fileprefix+"_N*.gsd"),key=lambda x:int(os.path.basename(x).split("_")[16].split('.')[0].replace('N',''))):
        #if(k>num_files-1):
            #break
        trajectory = gsd.open(filename,'rb') # read gsd file
        key=int(os.path.basename(filename).split("_")[16].split('.')[0].replace('N',''))
        print("Run no: ",key)
        print("No of frames: ",len(trajectory))
        """
        for i in range(len(trajectory)):
            if(k!=0 and i==0):
                continue
            system= trajectory[int(i)]
            traj.append(system) 
        k+=1
        """


