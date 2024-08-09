import os
import sys
import numpy as np
from gsd import hoomd as gsd
import argparse
import os.path
import glob

parser=argparse.ArgumentParser()
parser.add_argument("--fileprefix",type=str)

args = parser.parse_args()
locals().update(vars(args))

os.chdir(fileprefix)
for file in glob.glob("*.allruns.gsd"):
    outputfile=os.path.splitext(file)[0]+'_condensed.gsd'
    #print(outputfile)
    traj = gsd.open(name=outputfile,mode='wb')
    trajectory = gsd.open(file,'rb') #read gsd file
    print(len(trajectory))
    for i in range(0,len(trajectory),10):  #frame skipping
        system= trajectory[int(i)]
        traj.append(system)
    print(len(traj))

