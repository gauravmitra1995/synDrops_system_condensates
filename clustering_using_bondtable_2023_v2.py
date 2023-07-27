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

def remove_from_list(cluster_sizes,min_cluster_size):
    res=[]
    for c in cluster_sizes:
        if(c>=min_cluster_size):
            res.append(c)
    return res

def cluster_size_v_time(trajectory,minframe=0,frameinterval=1,verbose=False):
    
    num_frames=len(trajectory)
    maxframe=num_frames-1
    
    frames_list=np.arange(minframe,maxframe+1,frameinterval) # list of frames to analyze
    
    timestep_list=[]
    
    largest_clusters_with_time=[]
    #second_largest_clusters_with_time=[]
    numberofclusters_with_time=[]
    average_clustersize_with_time=[]

    print("Minimum cluster size: ",min_cluster_size)
    
    for frame in frames_list[::]:  #skipping no frames at all
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
        

        if(min_cluster_size==1):
            cluster_sizes.extend(singleton_clusters)  #add singletons to the list of cluster sizes too only if minimum cluster size chosen is 1 
        else:
            cluster_sizes_modified=remove_from_list(cluster_sizes,min_cluster_size)  #remove those cluster sizes which are less than the chosen minimum cluster size
            cluster_sizes=cluster_sizes_modified

        if len(cluster_sizes)>0:
            largest_cluster_size=cluster_sizes[0]
            
            #if(len(cluster_sizes)>1):
                #second_largest_cluster_size=cluster_sizes[1]
            #else:
                #second_largest_cluster_size=0
            
            average_cluster_size=np.round(np.mean(cluster_sizes),5)
        else:
            largest_cluster_size = 0 
            #second_largest_cluster_size=0
            average_cluster_size=0

        largest_clusters_with_time.append(largest_cluster_size)
        #second_largest_clusters_with_time.append(second_largest_cluster_size)
        average_clustersize_with_time.append(average_cluster_size)

        if verbose:
            print("*"*100)
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
    parser.add_argument("--trajectory_file",type=str,required=True)
    parser.add_argument("--min_cluster_size",type=int,required=True)
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
    
    largestclustersize_file='./largestclustersizevstime_data/'+os.path.splitext(os.path.basename(trajectory_file))[0]+'.largestclustersizevstime.data'
    averageclustersize_file='./averageclustersizevstime_data/'+os.path.splitext(os.path.basename(trajectory_file))[0]+'.averageclustersizevstime_minclustersize'+str(min_cluster_size)+'.data'

    max_cluster_size_data.dump(largestclustersize_file)
    average_cluster_size_data.dump(averageclustersize_file)
    
    
    

