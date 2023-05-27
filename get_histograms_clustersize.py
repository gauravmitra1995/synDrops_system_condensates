import os
import sys
import pickle
import numpy as np
import matplotlib.pyplot as plt
import glob
np.set_printoptions(suppress=True)
from gsd import hoomd as gsd
import random
from bond_analysis import get_bond_table
from collections import defaultdict
from clustering_using_bondtable_2023_v2 import *
from statistics import median


if __name__ == "__main__":
    import argparse
    parser=argparse.ArgumentParser()

    #parser.add_argument("--trajectory_file",type=str)
    parser.add_argument("--minframe",default=0,type=int)
    parser.add_argument("--frameinterval",default=1,type=int)
    parser.add_argument('--verbose',default=False,action='store_true')

    #Parameters which were varied for the conditions run here:
    parser.add_argument("--crowder_temperature",help="Temperature of crowders, relative to 1.0 (default: %(default)s)",default=1.0,type=float,required=True)
    parser.add_argument("--volume_fraction_ribosome",help="Fraction of volume taken up by 30nm ribosome (default: %(default)s)",default=0,type=float,required=True)
    parser.add_argument("--koff",help="Off rate for unbinding in time units (default: %(default)s)",default=0.001,type=float,required=True)
    
    #Parameters which were more or less fixed 
    parser.add_argument("--volume_fraction_polysome",help="Fraction of volume taken up by 100nm polysome (default: %(default)s)",default=0,type=int)
    parser.add_argument("--box_length",help="Side length of simulation box in nm (default: %(default)s)",default=860,type=float)
    parser.add_argument("--binding_distance",help="Maximum distance where binding can take place, in nm (default: %(default)s)",default=1.0,type=float)
    parser.add_argument("--dt",help="Simulation time step (default: %(default)s)",default=0.002,type=float)
    parser.add_argument("--nsteps",help="Number of steps of simulation to run (default: %(default)s)",default=100000000,type=int)

    parser.add_argument("--sphere_repulsion",help="Repulsion of spheres including crowders in units of kT, for soft potential (default: %(default)s)",default=500,type=float)
    parser.add_argument('--gamma_scale',help="Friction coefficient to multiply diameter to give diffusion coeff (default: %(default))",type=float,default=0.001)
    parser.add_argument("--number_rods",help="Number of rod proteins (default: %(default)s)",default=1170,type=int)
    parser.add_argument("--number_linkers",help="Number of rod proteins (default: %(default)s)",default=390,type=int)
    parser.add_argument("--number_gems",help="Number of GEM proteins (default: %(default)s)",default=20,type=int)

    args = parser.parse_args()
    locals().update(vars(args))

    checksteps=10
    kon=1/(dt*checksteps)

    if(koff==0.0):
        epsilon='infinite'
    else:
        if(kon!=0.0):
            epsilon=format(np.log(kon/koff),'.1f')
        else:
            epsilon=0.0

    max_cluster_sizes_allseeds=[]
    second_max_cluster_sizes_allseeds=[]
    average_cluster_sizes_allseeds=[]

    median_sizes_allseeds=[]

    time_allseeds=[]
    seedlist=[]

    #fig,ax=plt.subplots(figsize=(22,15),dpi=100)
    
    for filename in sorted(glob.glob(str(os.getcwd())+"/prod_v2.6_newdyn_2023/l"+str(box_length)+"_vfr"+str(volume_fraction_ribosome)+"_vfp"+str(volume_fraction_polysome)+"_Tc"+str(crowder_temperature)+"/gel_l"+str(box_length)+"_vfr"+str(volume_fraction_ribosome)+"_vfp"+str(volume_fraction_polysome)+"_nG"+str(number_gems)+"_nR"+str(number_rods)+"_nL"+str(number_linkers)+"_k00"+"_koff"+str(koff)+"_repuls"+str(sphere_repulsion)+"_bd"+str(binding_distance)+"_Tc"+str(crowder_temperature)+"_s*dt"+str(dt)+"_gs"+str(gamma_scale)+"_N"+str(nsteps)+".gsd"),key=lambda x:(int(((os.path.basename(x).split("_")[12]).split(".")[0]).replace('s','')))):
        seed=int(((os.path.basename(filename).split("_")[12]).split(".")[0]).replace('s',''))
        print("********************************************************************************************************************")
        print("Seed: ",seed)
        trajectory = gsd.open(filename,'rb') # read gsd file
        num_frames=len(trajectory)
        maxframe=num_frames-1
        frames_list=np.arange(minframe,maxframe+1,frameinterval) # list of frames to analyze
        print("Total number of frames: ",num_frames)
                
        if verbose:
            print('File:',trajectory_file)
            print('Total number of frames:',num_frames)

        max_cluster_size_data,second_max_cluster_size_data,average_cluster_size_data = cluster_size_v_time(trajectory,minframe,frameinterval,verbose)
        timesteps = max_cluster_size_data[:,0]
        step_time=7.5e-2
        time=timesteps*(step_time/1e6)
        time_allseeds.append(time)

        max_cluster_sizes=max_cluster_size_data[:,1]
        median_size=median(max_cluster_sizes)
        median_sizes_allseeds.append(median_size)

        second_max_cluster_sizes=second_max_cluster_size_data[:,1]

        average_cluster_sizes=average_cluster_size_data[:,1]

        #ax.plot(time,max_cluster_sizes,marker='o',linestyle=None,markersize=4.0,label='s = '+str(seed))

        max_cluster_sizes_allseeds.append(max_cluster_sizes)
        second_max_cluster_sizes_allseeds.append(second_max_cluster_sizes)
        average_cluster_sizes_allseeds.append(average_cluster_sizes)

        seedlist.append(seed)

    median_cluster_size = [median(frame) for frame in zip(*max_cluster_sizes_allseeds)]
    mean_cluster_size = [float(sum(frame)/len(frame)) for frame in zip(*max_cluster_sizes_allseeds)]
    shortest_trajectory = min(time_allseeds, key=len)

    median_of_medians_framebyframe = median(median_cluster_size)
    median_of_medians_seedbyseed = median(median_sizes_allseeds)

    #print("Median of medians frame by frame: ",median_of_medians_framebyframe)
    #print("Median of medians seed by seed: ",median_of_medians_seedbyseed)
    """
    ax.plot(shortest_trajectory,mean_cluster_size,linestyle='dashed',linewidth=5.0,color='black',label='Mean')
    ax.plot(shortest_trajectory,median_cluster_size,linestyle='dotted',linewidth=5.0,color='black',label='Median')

    ax.set_ylabel('Largest cluster size',fontsize=60)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(5)
    ax.set_xlabel('Time (in sec)',fontsize=60)
    plt.xticks(fontsize=50)
    plt.yticks(fontsize=50)
    plt.title(r'$\phi_{ribo} = $'+str(volume_fraction_ribosome)+' ; '+r'$T_{c} = $'+str(crowder_temperature)+' ; '+r'$\varepsilon = $'+str(epsilon),fontsize=60)
    ax.legend(loc='upper left',ncol=1,prop={'size': 40})
    fig.tight_layout()
    plt.savefig('final_figures/largestclustersize_vs_time/volfracribo'+str(volume_fraction_ribosome)+'_Tc'+str(crowder_temperature)+'_eps'+str(epsilon)+'_largestclustersizevstime.png',bbox_inches='tight')
    plt.close()
    """

    fig,ax=plt.subplots(figsize=(22,15),dpi=100)

    k=0
    for a in average_cluster_sizes_allseeds:
        length=a.shape[0]
        seed=seedlist[k]
        time=time_allseeds[k]
        print(time[0],a[0])
        ax.plot(time[:],a[:],marker='o',linewidth=5.0,markersize=10.0,label='s = '+str(seed))
        k+=1
    ax.set_ylabel('Average cluster size',fontsize=60)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(5)
    ax.set_xlabel('Time (in sec)',fontsize=60)
    plt.xticks(fontsize=50)
    plt.yticks(fontsize=50)
    plt.title(r'$\phi_{ribo} = $'+str(volume_fraction_ribosome)+' ; '+r'$T_{c} = $'+str(crowder_temperature)+' ; '+r'$\varepsilon = $'+str(epsilon),fontsize=60)
    plt.yscale('log')
    plt.xscale('log')
    ax.legend(loc='upper left',ncol=1,prop={'size': 40})
    fig.tight_layout()
    plt.xlim(0.1,10)
    #plt.ylim(,None)
    plt.savefig('final_figures/averageclustersize_vs_time/volfracribo'+str(volume_fraction_ribosome)+'_Tc'+str(crowder_temperature)+'_eps'+str(epsilon)+'_averageclustersizevstime.png',bbox_inches='tight')
    plt.close()


    """   
    fig,ax=plt.subplots(figsize=(22,15),dpi=100)

    k=0
    for m in max_cluster_sizes_allseeds:
        length=m.shape[0]
        seed=seedlist[k]
        ax.hist(m[int(8*length/10.):],bins=50,label='s = '+str(seed),width=5,density=True,alpha=0.5)
        k+=1

    ax.legend(loc='upper right',ncol=2,prop={'size': 40})
    ax.set_xlabel('Largest cluster size',fontsize=60)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(5)
    ax.set_ylabel('Probability',fontsize=60)
    plt.xticks(fontsize=50)
    plt.yticks(fontsize=50)
    ax.set_xlim(0,1500)
    plt.title(r'$\phi_{ribo} = $'+str(volume_fraction_ribosome)+' ; '+r'$T_{c} = $'+str(crowder_temperature)+' ; '+r'$\varepsilon = $'+str(epsilon),fontsize=60)
    fig.tight_layout()
    plt.savefig('final_figures/histograms/volfracribo'+str(volume_fraction_ribosome)+'_Tc'+str(crowder_temperature)+'_eps'+str(epsilon)+'_largestclustersizehistogram.png',bbox_inches='tight')
    plt.close()    
    """



    """
    fig,ax=plt.subplots(figsize=(22,15),dpi=100)

    x=0
    for sm in second_max_cluster_sizes_allseeds:
        length=sm.shape[0]
        seed=seedlist[x]
        ax.hist(sm[int(8*length/10.):],bins=50,label='s = '+str(seed),width=5,density=True,alpha=0.5)
        x+=1

    ax.legend(loc='upper right',ncol=2,prop={'size': 40})
    ax.set_xlabel('Second largest cluster size',fontsize=60)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(5)
    ax.set_ylabel('Probability',fontsize=60)
    plt.xticks(fontsize=50)
    plt.yticks(fontsize=50)
    ax.set_xlim(0,1500)
    plt.title(r'$\phi_{ribo} = $'+str(volume_fraction_ribosome)+' ; '+r'$T_{c} = $'+str(crowder_temperature)+' ; '+r'$\varepsilon = $'+str(epsilon),fontsize=60)
    fig.tight_layout()
    plt.savefig('final_figures/histograms/volfracribo'+str(volume_fraction_ribosome)+'_Tc'+str(crowder_temperature)+'_eps'+str(epsilon)+'_secondlargestclustersizehistogram.png',bbox_inches='tight')
    plt.close()
    """
        
        
    

