import os
import sys
import pickle
import numpy as np
import matplotlib.pyplot as plt
np.set_printoptions(suppress=True)
from gsd import hoomd as gsd
import random
from bond_analysis import get_bond_table
from collections import defaultdict

class UnionFind:
    def __init__(self):
        self.parent = {}

    def find(self, x):
        if x not in self.parent:
            self.parent[x] = x
        elif self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, x, y):
        px, py = self.find(x), self.find(y)
        if px != py:
            self.parent[px] = py

def get_clustering_fast(cluster_bondtable_new):
    uf = UnionFind()
    for pair in cluster_bondtable_new:
        uf.union(pair[0], pair[1])

    clusters = defaultdict(set)
    for pair in cluster_bondtable_new:
        root = uf.find(pair[0])
        clusters[root].update(pair)

    output = [sorted(list(cluster)) for cluster in clusters.values()]

    return output

def flat_list(alist): # [[1,2],[3,4]] --> [1,2,3,4]
    flatlist=list()
    for sublist in alist:
        for i in sublist:
            flatlist.append(i)
    return(flatlist)

def find_singletons(combined_list,all_bodies):
    setbigger = set(all_bodies)
    setsmaller = set(combined_list)
    missing_elements = setbigger - setsmaller
    return list(missing_elements)


def cluster_size_v_time(trajectory,minframe=0,frameinterval=1,verbose=False):
    
    num_frames=len(trajectory)
    maxframe=num_frames-1
    
    frames_list=np.arange(minframe,maxframe+1,frameinterval) # list of frames to analyze
    
    timestep_list=[]
    
    largest_clusters_with_time=[]
    second_largest_clusters_with_time=[]
    numberofclusters_with_time=[]
    average_clustersize_with_time=[]
    
    for frame in frames_list[::]:
        system=trajectory[int(frame)]

        cluster_bondtable,bodies = get_bond_table(trajectory,frame_id=int(frame))
        all_bodies=list(np.arange(0,np.max(bodies)+1))

        clustering = get_clustering_fast(cluster_bondtable)
        combined_list=sorted(flat_list(clustering))

        singletons=find_singletons(combined_list,all_bodies)
        number_singletons=len(singletons)
        singleton_clusters=list(np.ones(number_singletons).astype(int))
            
        numberofclusters_with_time.append(len(clustering))
        timestep_list.append(system.configuration.step)
        cluster_sizes=sorted([len(c) for c in clustering],reverse=True)
        cluster_sizes.extend(singleton_clusters)

        if len(cluster_sizes)>0:
            largest_cluster_size=cluster_sizes[0]
            if(len(cluster_sizes)>1):
                second_largest_cluster_size=cluster_sizes[1]
            else:
                second_largest_cluster_size=0
            average_cluster_size=np.round(np.mean(cluster_sizes),5)
        else:
            largest_cluster_size = 0 
            second_largest_cluster_size=0
            average_cluster_size=0

        largest_clusters_with_time.append(largest_cluster_size)
        second_largest_clusters_with_time.append(second_largest_cluster_size)
        average_clustersize_with_time.append(average_cluster_size)

        if verbose:
            print("*"*80)
            print("Frame:",frame)
            print("No of clusters:",len(cluster_sizes))
            print(cluster_sizes)
            print("Total particles:",np.sum(cluster_sizes))
            print("Largest cluster size:",largest_cluster_size)
            print("Average cluster size:",average_cluster_size)
            #print("Second largest cluster size:",second_largest_cluster_size)
        
    
    largest_clusters_with_time = np.array(largest_clusters_with_time).reshape(-1,1)
    #second_largest_clusters_with_time = np.array(second_largest_clusters_with_time).reshape(-1,1)
    average_clustersize_with_time = np.array(average_clustersize_with_time).reshape(-1,1)


    timestep_list = np.array(timestep_list).reshape(-1,1)

    max_cluster_size_data = np.concatenate((timestep_list,largest_clusters_with_time),axis=1)
    #second_max_cluster_size_data = np.concatenate((timestep_list,second_largest_clusters_with_time),axis=1)
    average_cluster_size_data = np.concatenate((timestep_list,average_clustersize_with_time),axis=1)

    return max_cluster_size_data,average_cluster_size_data
    
    

if __name__ == "__main__":
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument("--trajectory_file",type=str)
    parser.add_argument("--minframe",default=0,type=int)
    parser.add_argument("--frameinterval",default=1,type=int)
    parser.add_argument('--verbose',default=False,action='store_true')

    args = parser.parse_args()
    locals().update(vars(args))
    
    trajectory = gsd.open(trajectory_file,'rb') # read gsd file
    if verbose:
        print('File:',trajectory_file)
        print('Total number of frames:',len(trajectory))

    max_cluster_size_data,average_cluster_size_data  = cluster_size_v_time(trajectory,minframe,frameinterval,verbose)
    
    max_cluster_sizes=max_cluster_size_data[:,1]
    average_cluster_sizes=average_cluster_size_data[:,1]

    timesteps=max_cluster_size_data[:,0]
    step_time=7.5e-2
    time=timesteps*(step_time/1e6)


    """   
    fig,ax=plt.subplots(figsize=(22,15),dpi=100)

    #Linear fitting the log scale graph (check power law index)

    nucleation_time_guess = 0.5

    coef = np.polyfit(np.log10(time),np.log10(average_cluster_sizes),1)
    poly1d_fn = np.poly1d(coef) # poly1d_fn is now a function which takes in x and returns an estimate for y
    cluster_time_fit = [pow(10,i) for i in poly1d_fn(np.log10(time))]
    ax.plot(time,average_cluster_sizes,linewidth=5.0,color='blue')
    ax.plot(time,cluster_time_fit,'--r',label=f'Total: Exponent $\\alpha$={round(coef[0], 2)}')

    first_select = np.where(time < nucleation_time_guess)[0]

    coef = np.polyfit(np.log10(time[first_select]),np.log10(average_cluster_sizes[first_select]),1)
    poly1d_fn = np.poly1d(coef)
    cluster_time_fit = [pow(10,i) for i in poly1d_fn(np.log10(time))]
    ax.plot(time,cluster_time_fit,'--g',label=f'First: Exponent $\\alpha$={round(coef[0], 2)}')

    #Make the plots

    #ax.plot(time,max_cluster_sizes,marker='o',linewidth=5.0,markersize=10.0,color='blue')

    ax.set_ylabel('Average Cluster size',fontsize=60)
    #ax.set_ylabel('Largest Cluster size',fontsize=60)
    ax.set_xlabel('Time (in sec)',fontsize=60)

    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(5)
    ax.legend(loc='upper left',ncol=1,prop={'size': 40})
    plt.xticks(fontsize=50)
    plt.yticks(fontsize=50)
    plt.yscale('log')
    plt.xscale('log')
    fig.tight_layout()
    plt.show()
    plt.close()
    """
    
    
    


        


