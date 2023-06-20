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


if __name__ == "__main__":
    import argparse
    parser=argparse.ArgumentParser()

    #Parameters which were varied for the conditions run here:
    parser.add_argument("--crowder_temperature",help="Temperature of crowders, relative to 1.0 (default: %(default)s)",default=1.0,type=float,required=True)
    parser.add_argument("--volume_fraction_ribosome",help="Fraction of volume taken up by 30nm ribosome (default: %(default)s)",default=0,type=float,required=True)
    parser.add_argument("--koff",help="Off rate for unbinding in time units (default: %(default)s)",default=0.001,type=float,required=True)
    parser.add_argument("--min_cluster_size",help="Minimum cluster size to be considered",default=1,type=int,required=True)
    
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

    if(kon==0.0):
        epsilon=0.0
    else:
        if(koff!=0.0):
            epsilon=format(np.log(kon/koff),'.1f')
        else:
            epsilon='infinite'


    print("/"*100)
    print("Epsilon under consideration: ",epsilon)
    print("Volume fraction of ribosome: ",volume_fraction_ribosome)
    print("Crowder temperature: ",crowder_temperature)
    print("/"*100)

    '''        
    #LARGEST CLUSTER SIZE VS TIME ANALYSIS
    #For largest cluster size vs time, combine all seeds and put them individually on same plot

    time_allseeds=[]
    seedlist=[]
    max_cluster_sizes_allseeds=[]
    #median_sizes_allseeds=[]

    fig,ax=plt.subplots(figsize=(20,15),dpi=100)


    cutoff_frames=2400
    for filename in sorted(glob.glob("./largestclustersizevstime_data/gel_l"+str(box_length)+"_vfr"+str(volume_fraction_ribosome)+"_vfp"+str(volume_fraction_polysome)+"_nG"+str(number_gems)+"_nR"+str(number_rods)+"_nL"+str(number_linkers)+"_k0"+str(koff0)+"_koff"+str(koff)+"_repuls"+str(sphere_repulsion)+"_bd"+str(binding_distance)+"_Tc"+str(crowder_temperature)+"_s*dt"+str(dt)+"_gs"+str(gamma_scale)+"_combined.largestclustersizevstime.data"),key=lambda x:(int(((os.path.basename(x).split("_")[12]).split(".")[0]).replace('s','')))):
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
        print(max_cluster_sizes)
        print(max_cluster_sizes.shape[0])

        #median_size=median(max_cluster_sizes)

        ax.plot(time,max_cluster_sizes,marker=None,linestyle='-',linewidth=3.0,label='s = '+str(seed))

        max_cluster_sizes_allseeds.append(max_cluster_sizes)
        #median_sizes_allseeds.append(median_size)

        seedlist.append(seed)

    """
    median_cluster_size = [median(frame) for frame in zip(*max_cluster_sizes_allseeds)]
    mean_cluster_size = [float(sum(frame)/len(frame)) for frame in zip(*max_cluster_sizes_allseeds)]
    shortest_trajectory = min(time_allseeds, key=len)

    median_of_medians_framebyframe = median(median_cluster_size)
    median_of_medians_seedbyseed = median(median_sizes_allseeds)

    #print("Median of medians frame by frame: ",median_of_medians_framebyframe)
    #print("Median of medians seed by seed: ",median_of_medians_seedbyseed)
    """

    #ax.plot(shortest_trajectory,mean_cluster_size,linestyle='dashed',linewidth=5.0,color='black',label='Mean')
    #ax.plot(shortest_trajectory,median_cluster_size,linestyle='dotted',linewidth=5.0,color='black',label='Median')


    
    ax.set_ylabel('Largest cluster size',fontsize=50)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(5)
    ax.set_xlabel('Time (in sec)',fontsize=50)
    ax.set_ylim(0,1600)
    plt.xticks(np.arange(0,26.0,2.0),fontsize=40)
    plt.yticks(np.arange(0,1600,100),fontsize=40)
    plt.grid()
    plt.title(r'$\phi_{ribosome} = $'+str(volume_fraction_ribosome)+' ; '+r'$T_{c} = $'+str(crowder_temperature)+' ; '+r'$\varepsilon = $'+str(epsilon),fontsize=50)
    ax.legend(loc='upper left',ncol=1,prop={'size': 30})
    fig.tight_layout()
    plt.savefig('final_figures/largestclustersize_vs_time/volfracribo'+str(volume_fraction_ribosome)+'_Tc'+str(crowder_temperature)+'_eps'+str(epsilon)+'_largestclustersizevstime.svg',bbox_inches='tight')
    plt.close()

    '''
    

    
    #AVERAGE CLUSTER SIZE VS TIME ANALYSIS
    #For largest cluster size vs time, combine all seeds by averaging them and put the average on the plot with error bars indicating standard deviation.

    time_allseeds=[]
    seedlist=[]
    average_cluster_sizes_allseeds=[]
    
    fig,ax=plt.subplots(figsize=(20,15),dpi=100)

    cutoff_frames=2400

    for filename in sorted(glob.glob("./averageclustersizevstime_data/gel_l"+str(box_length)+"_vfr"+str(volume_fraction_ribosome)+"_vfp"+str(volume_fraction_polysome)+"_nG"+str(number_gems)+"_nR"+str(number_rods)+"_nL"+str(number_linkers)+"_k0"+str(koff0)+"_koff"+str(koff)+"_repuls"+str(sphere_repulsion)+"_bd"+str(binding_distance)+"_Tc"+str(crowder_temperature)+"_s*dt"+str(dt)+"_gs"+str(gamma_scale)+"_combined.averageclustersizevstime_minclustersize"+str(min_cluster_size)+".data"),key=lambda x:(int(((os.path.basename(x).split("_")[12]).split(".")[0]).replace('s','')))):
        seed=int(((os.path.basename(filename).split("_")[12]).split(".")[0]).replace('s',''))
        #print("*"*100)
        #print("Seed: ",seed)

        average_cluster_size_data=np.load(filename,allow_pickle=True)

        timesteps = average_cluster_size_data[:,0]
        step_time=7.5e-2
        time=timesteps*(step_time/1e6)

        print("No of frames in combined gsd before truncation: ",time.shape[0])

        time=time[:cutoff_frames]

        print("No of frames in combined gsd after truncation: ",time.shape[0])
        time_allseeds.append(time)

        average_cluster_sizes=average_cluster_size_data[:,1].astype(float)[:cutoff_frames]
        print(average_cluster_sizes)
        print(average_cluster_sizes.shape[0])

        average_cluster_sizes_allseeds.append(average_cluster_sizes)

        seedlist.append(seed)
  
    print("Length of the data for the 5 seeds: ")
    print([len(i) for i in average_cluster_sizes_allseeds])

    #for averaging make sure the seeds are of the same length (#frames)
    min_length=min([len(i) for i in average_cluster_sizes_allseeds])

    timesteps=average_cluster_size_data[:,0][:min_length]
    step_time=7.5e-2
    time=timesteps*(step_time/1e6)

    average_cluster_sizes_allseeds_allreshaped=[]
    for a in average_cluster_sizes_allseeds:
        a=a[:min_length]
        average_cluster_sizes_allseeds_allreshaped.append(a)

    average_cluster_sizes_allseeds_allreshaped=np.array(average_cluster_sizes_allseeds_allreshaped,dtype=object)

    average_cluster_sizes_meanofallseeds=np.mean(average_cluster_sizes_allseeds_allreshaped.astype(float),axis=0)
    average_cluster_sizes_stddevofallseeds=np.std(average_cluster_sizes_allseeds_allreshaped.astype(float),axis=0)

    #Linear fitting the log scale graph (check power law index)

    nucleation_time_guess = 0.5

    coef = np.polyfit(np.log10(time),np.log10(average_cluster_sizes_meanofallseeds),1)
    poly1d_fn = np.poly1d(coef) # poly1d_fn is now a function which takes in x and returns an estimate for y
    cluster_time_fit = [pow(10,i) for i in poly1d_fn(np.log10(time))]
    ax.errorbar(time[::10],average_cluster_sizes_meanofallseeds[::10],yerr=average_cluster_sizes_stddevofallseeds[::10],marker=None,linestyle='-',color='blue',linewidth=5.0,ecolor='k',elinewidth=1.0,capsize=5,capthick=1.0)
    ax.plot(time,cluster_time_fit,'--r',label=f'Total: Exponent $\\alpha$={round(coef[0], 2)}')

    first_select = np.where(time < nucleation_time_guess)[0]

    coef = np.polyfit(np.log10(time[first_select]),np.log10(average_cluster_sizes_meanofallseeds[first_select]),1)
    poly1d_fn = np.poly1d(coef)
    cluster_time_fit = [pow(10,i) for i in poly1d_fn(np.log10(time))]
    ax.plot(time,cluster_time_fit,'--g',label=f'First: Exponent $\\alpha$={round(coef[0], 2)}')

    #Make the plots

    ax.set_ylabel('Average Cluster size',fontsize=50)
    ax.set_xlabel('Time (in sec)',fontsize=50)
    ax.set_ylim(1,100)

    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(5)
    ax.legend(loc='upper left',ncol=1,prop={'size': 30})
    plt.xticks(fontsize=40)
    plt.yticks(fontsize=40)
    plt.yscale('log')
    plt.xscale('log')
    plt.title(r'$\phi_{ribosome} = $'+str(volume_fraction_ribosome)+' ; '+r'$T_{c} = $'+str(crowder_temperature)+' ; '+r'$\varepsilon = $'+str(epsilon)+' ; '+'Min clus size = '+str(min_cluster_size),fontsize=40)
    fig.tight_layout()
    plt.savefig('final_figures/averageclustersize_vs_time/volfracribo'+str(volume_fraction_ribosome)+'_Tc'+str(crowder_temperature)+'_eps'+str(epsilon)+'_averageclustersizevstime_minclustersize'+str(min_cluster_size)+'.svg',bbox_inches='tight')
    plt.close()
    

    
    




