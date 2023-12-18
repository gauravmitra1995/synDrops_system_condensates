import os
import sys
import pickle
import numpy as np
import matplotlib.pyplot as plt
import glob
np.set_printoptions(suppress=True)
from gsd import hoomd as gsd
import random
from collections import defaultdict
from statistics import median
from functions_cluster_analysis import *
from matplotlib.lines import Line2D
import seaborn as sns
sns.set_context("talk", font_scale=1.0)



if __name__ == "__main__":
    import argparse
    parser=argparse.ArgumentParser()

    #Parameters which were varied for the conditions run here:
    parser.add_argument("--crowder_temperature",help="Temperature of crowders, relative to 1.0 (default: %(default)s)",default=1.0,type=float)
    #parser.add_argument("--volume_fraction_ribosome",help="Fraction of volume taken up by 30nm ribosome (default: %(default)s)",default=0.0,type=float,required=True)
    #parser.add_argument("--koff",help="Off rate for unbinding in time units (default: %(default)s)",default=0.001,type=float,required=True)
    
    #Parameters which were more or less fixed 
    parser.add_argument("--volume_fraction_polysome",help="Fraction of volume taken up by 100nm polysome (default: %(default)s)",default=0,type=int)
    parser.add_argument("--box_length",help="Side length of simulation box in nm (default: %(default)s)",default=860,type=float)
    parser.add_argument("--binding_distance",help="Maximum distance where binding can take place, in nm (default: %(default)s)",default=1.0,type=float)
    parser.add_argument("--dt",help="Simulation time step (default: %(default)s)",default=0.002,type=float)
    parser.add_argument("--koff0",default=0,type=int)

    parser.add_argument("--sphere_repulsion",help="Repulsion of spheres including crowders in units of kT, for soft potential (default: %(default)s)",default=500,type=float)
    parser.add_argument('--gamma_scale',help="Friction coefficient to multiply diameter to give diffusion coeff (default: %(default))",type=float,default=0.001)
    parser.add_argument("--number_rods",help="Number of rod proteins (default: %(default)s)",default=1170,type=int)
    parser.add_argument("--number_linkers",help="Number of rod proteins (default: %(default)s)",default=390,type=int)
    parser.add_argument("--number_gems",help="Number of GEM proteins (default: %(default)s)",default=20,type=int)

    args = parser.parse_args()
    locals().update(vars(args))

    checksteps=10
    kon=1/(dt*checksteps)
 
    #pretty_print=lambda x: np.format_float_positional(x, trim="-")
    #koff=pretty_print(koff)
 
    print("/"*100)
    #print("Epsilon under consideration: ",epsilon)
    #print("Volume fraction of ribosome: ",volume_fraction_ribosome)
    print("Crowder temperature: ",crowder_temperature)
    #print("koff: ",koff)
    print("/"*100)

         
    #PLOT MEDIAN LARGEST CLUSTER SIZE VS EPSILON

    fig,ax=plt.subplots(figsize=(7,5.5),dpi=300)

    epsilons=[]
    medianlargestclustersizes_dict={}
    vfrlist=[]

    for filename in sorted(glob.glob("./../largestclustersizevstime_data/epsilonvariation/volfracribo*_Tc"+str(crowder_temperature)+"_eps*_medianlargestclustersize.data"),key=lambda x:(float((os.path.basename(x).split("_")[0]).replace('volfracribo','')),float((os.path.basename(x).split("_")[2]).replace('eps','')))):
        
        data=np.load(filename,allow_pickle=True)
        vfr=float((os.path.basename(filename).split("_")[0]).replace('volfracribo',''))
        epsilon=float(data[0])
        label=(vfr,epsilon)
        largestclustersize=data[1]
        #print(label)   
        if(epsilon not in epsilons):
            epsilons.append(float(epsilon))
        if vfr not in medianlargestclustersizes_dict:
            medianlargestclustersizes_dict[vfr] = [largestclustersize]
        else:
            medianlargestclustersizes_dict[vfr].append(largestclustersize)
  
    print(medianlargestclustersizes_dict)
    differences = [abs(x - y) for x, y in zip(medianlargestclustersizes_dict[0.0], medianlargestclustersizes_dict[0.3])]
    print(differences)
    print(epsilons)
    
    ax.plot(epsilons[:-1],differences[:-1],marker='o',linestyle='-',markersize=15.0,linewidth=5.0,color='k')
    ax.set_ylabel('Difference in largest cluster size')
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(2)
    ax.set_xlabel('Binding affinity ('+r'$\varepsilon$'+') in '+r'$k_B T$')
    ax.set_ylim(0,700)
    plt.xticks(np.arange(int(np.min(epsilons[:-1])),int(np.max(epsilons[:-1]))+1,2))
    plt.yticks(np.arange(0,800,100))
    plt.axvline(epsilons[differences.index(np.max(differences))],linestyle='--',linewidth=3.0,color='g')
    #plt.title(r'$\phi_{ribosome} = $'+str(volume_fraction_ribosome)+' ; '+r'$T_{c} = $'+str(crowder_temperature))
    #plt.grid(alpha=0.4)
    fig.tight_layout()
    plt.savefig('./../final_figures/largestclustersize_vs_time/epsilonvariation/Tc'+str(crowder_temperature)+'_differenceinlargestclustersizevsepsilon_between2vfr.pdf',bbox_inches='tight')
    plt.close()
    
    
   
