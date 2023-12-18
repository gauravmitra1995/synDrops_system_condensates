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


if __name__ == "__main__":
    import argparse
    parser=argparse.ArgumentParser()

    #Parameters which were varied for the conditions run here:
    parser.add_argument("--crowder_temperature",help="Temperature of crowders, relative to 1.0 (default: %(default)s)",default=1.0,type=float,required=True)
    parser.add_argument("--volume_fraction_ribosome",help="Fraction of volume taken up by 30nm ribosome (default: %(default)s)",default=0,type=float,required=True)
    parser.add_argument("--koff",help="Off rate for unbinding in time units (default: %(default)s)",default=0.001,type=float,required=True)
    parser.add_argument("--number_rods",help="Number of rod proteins (default: %(default)s)",default=1170,type=int,required=True)
    parser.add_argument("--number_linkers",help="Number of rod proteins (default: %(default)s)",default=390,type=int,required=True)
    parser.add_argument("--min_cluster_size",help="Minimum cluster size to be considered",default=1,type=int,required=True)
    parser.add_argument("--cutoff_frames",help="Cut off value for the frames upto which analysis has to be done",default=2400,type=int,required=True)
    
    #Parameters which were more or less fixed 
    parser.add_argument("--volume_fraction_polysome",help="Fraction of volume taken up by 100nm polysome (default: %(default)s)",default=0,type=int)
    parser.add_argument("--box_length",help="Side length of simulation box in nm (default: %(default)s)",default=860,type=float)
    parser.add_argument("--binding_distance",help="Maximum distance where binding can take place, in nm (default: %(default)s)",default=1.0,type=float)
    parser.add_argument("--dt",help="Simulation time step (default: %(default)s)",default=0.002,type=float)
    parser.add_argument("--koff0",default=0,type=int)

    parser.add_argument("--sphere_repulsion",help="Repulsion of spheres including crowders in units of kT, for soft potential (default: %(default)s)",default=500,type=float)
    parser.add_argument('--gamma_scale',help="Friction coefficient to multiply diameter to give diffusion coeff (default: %(default))",type=float,default=0.001)
    parser.add_argument("--number_gems",help="Number of GEM proteins (default: %(default)s)",default=20,type=int)

    args = parser.parse_args()
    locals().update(vars(args))

    checksteps=10
    kon=1/(dt*checksteps)

    if(kon==0.0):
        epsilon=0.0
    else:
        if(koff!=0.0):
            epsilon=format(np.log(kon/koff),'.1f')
        else:
            epsilon='infinite'

    
    pretty_print=lambda x: np.format_float_positional(x, trim="-")
    koff=pretty_print(koff)

    print("/"*100)
    print("Epsilon under consideration: ",epsilon)
    print("Volume fraction of ribosome: ",volume_fraction_ribosome)
    print("Crowder temperature: ",crowder_temperature)
    print("koff: ",koff)
    print("No of rods: ",number_rods)
    print("No of linkers: ",number_linkers)
    print("/"*100)

         
    #LARGEST CLUSTER SIZE ANALYSIS

    time_allseeds=[]
    seedlist=[]
    max_cluster_sizes_allseeds=[]

    max_cluster_size_final=[]
    max_cluster_size_9seconds=[]

    fig,ax=plt.subplots(figsize=(20,15),dpi=100)

    cutoff_frames=2400  

    for filename in sorted(glob.glob("./../largestclustersizevstime_data/gel_l"+str(box_length)+"_vfr"+str(volume_fraction_ribosome)+"_vfp"+str(volume_fraction_polysome)+"_nG"+str(number_gems)+"_nR"+str(number_rods)+"_nL"+str(number_linkers)+"_k0"+str(koff0)+"_koff"+str(koff)+"_repuls"+str(sphere_repulsion)+"_bd"+str(binding_distance)+"_Tc"+str(crowder_temperature)+"_s*dt"+str(dt)+"_gs"+str(gamma_scale)+".allruns.largestclustersizevstime.data"),key=lambda x:(int(((os.path.basename(x).split("_")[12]).split(".")[0]).replace('s','')))):
        seed=int(((os.path.basename(filename).split("_")[12]).split(".")[0]).replace('s',''))
        print("*"*100)
        print("Seed: ",seed)
        max_cluster_size_data=np.load(filename,allow_pickle=True)
        timesteps = max_cluster_size_data[:,0]
        step_time=7.5e-2
        time=timesteps*(step_time/1e6)
        
        print("No of frames in combined gsd before truncation: ",time.shape[0])
        time=time[:cutoff_frames]

        print("No of frames in combined gsd after truncation: ",time.shape[0])
        time_allseeds.append(time)

        max_cluster_sizes=max_cluster_size_data[:,1].astype(int)[:cutoff_frames]

        choice_of_time=9.0
        abs_diff = np.abs(time - choice_of_time)

        # Find the index of the minimum absolute difference
        index_of_closest_value = np.argmin(abs_diff)
        
        print(time[index_of_closest_value],max_cluster_sizes[index_of_closest_value])
        print(time[-1],max_cluster_sizes[-1])

        ax.plot(time,max_cluster_sizes,marker=None,linestyle='-',linewidth=3.0,label='s = '+str(seed))
  
        max_cluster_sizes_allseeds.append(max_cluster_sizes)
        max_cluster_size_final.append(max_cluster_sizes[-1])
        max_cluster_size_9seconds.append(max_cluster_sizes[index_of_closest_value])

        seedlist.append(seed) 
    

    median_max_cluster_size_final=median(max_cluster_size_final)
    median_max_cluster_size_9seconds=median(max_cluster_size_9seconds)

    #mean_max_cluster_size_final=np.mean(max_cluster_size_final,axis=0)
    #stddev_max_cluster_size_final=np.std(max_cluster_size_final,axis=0)

    print(number_rods,number_linkers,median_max_cluster_size_final,median_max_cluster_size_9seconds)
        
    output_file='./largestclustersizevstime_data/hexamerconcvariation/volfracribo'+str(volume_fraction_ribosome)+'_nR'+str(number_rods)+'_nL'+str(number_linkers)+'_Tc'+str(crowder_temperature)+'_eps'+str(epsilon)+'_medianlargestclustersize.data'
    d=np.array([number_rods,number_linkers,median_max_cluster_size_final,median_max_cluster_size_9seconds],dtype=object)
    d.dump(output_file)
    
  
    ax.set_ylabel('Largest cluster size',fontsize=50)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(5)
    ax.set_xlabel('Time (in sec)',fontsize=50)
    #ax.set_ylim(0,(number_rods+number_linkers)+100)
    plt.xticks(np.arange(0,26.0,2.0),fontsize=40)
    plt.yticks(fontsize=40)
    plt.grid(alpha=0.6)
    plt.title(r'$\phi_{ribosome} = $'+str(volume_fraction_ribosome)+' ; '+r'$N_{l} = $'+str(number_linkers)+' ; '+r'$T_{c} = $'+str(crowder_temperature)+' ; '+r'$\varepsilon = $'+str(epsilon),fontsize=50)
    ax.legend(loc='upper left',ncol=1,prop={'size': 30})
    fig.tight_layout()
    plt.savefig('./../final_figures/largestclustersize_vs_time/hexamerconcvariation/volfracribo'+str(volume_fraction_ribosome)+'_nR'+str(number_rods)+'_nL'+str(number_linkers)+'_Tc'+str(crowder_temperature)+'_eps'+str(epsilon)+'_largestclustersizevstime.svg',bbox_inches='tight')
    plt.close()
    
   
