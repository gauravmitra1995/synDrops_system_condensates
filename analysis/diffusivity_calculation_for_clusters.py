import os
import sys
import pickle
import numpy as np
import matplotlib.pyplot as plt
np.set_printoptions(suppress=True)
from gsd import hoomd as gsd
import random
from bond_analysis import get_bond_table
from bond_analysis import bonds_analysis
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

def group_elements(elements):
    grouped_elements = []  # The list to store the grouped elements
    current_group = []  # The current group being formed

    for element in elements:
        if not current_group or element == current_group[-1]:
            # If the current group is empty or the element is the same as the last element in the group
            current_group.append(element)  # Add the element to the current group
        else:
            # If the element is different from the last element in the group
            grouped_elements.append(current_group)  # Add the current group to the grouped_elements list
            current_group = [element]  # Start a new group with the current element

    grouped_elements.append(current_group)  # Add the last group to the grouped_elements list

    return grouped_elements


def group_elements_basedonalist(elements1,elements2):
    grouped_elements = []
    current_group = []
    last_element = None

    for element1,element2  in zip(elements1,elements2):
        if element1 != last_element:
            if current_group:
                grouped_elements.append(current_group)
            current_group = []
        current_group.append(element2)
        last_element = element1

    if current_group:
        grouped_elements.append(current_group)

    return grouped_elements


def get_trajectoryinfo_vs_time(trajectory,cutoff_frames,minframe,frameinterval,verbose=False,step_time=7.5e-2):
    
    num_frames=len(trajectory)
    maxframe=num_frames-1
    
    frames_list=np.arange(minframe,maxframe+1,frameinterval) # list of frames to analyze
    timestep_list=[]
    time_list=[]

    box_length = None
    types = trajectory[0].particles.typeid
    N_A = int((types==0).sum())
    N_B = int((types==1).sum())
    n_binders = N_A+N_B
    assert N_A > N_B, "Warning, these may not be the right labels for A (n=%i) and B (n=%i)"%(N_A,N_B)

    prev_positions = np.zeros((n_binders,3))
    deltaT_seconds = step_time/1e6

    print("Delta time in seconds: ",deltaT_seconds)
    print("Minimum cluster size: ",min_cluster_size)
    print("Cut off frames for analysis: ",cutoff_frames)
    print("Frame interval: ",frameinterval)

    frames_list=[f for f in frames_list if f<(cutoff_frames+1)]

    #print(np.array(frames_list))
    print("No of frames for analysis: ",len(frames_list))

    mean_diffusivities_allframes=[]
    std_diffusivities_allframes=[]
    cluster_sizes_allframes=[]

    nucleation_time_guess = 0.5
    
    largestclustersizes=[]
    

    for frame in frames_list:   #doing analysis only upto cutoff frames
        system=trajectory[int(frame)]
        step=system.configuration.step

        box_size_i = trajectory[int(frame)].configuration.box[:3]
        if box_length is None:
            box_length = box_size_i[0]
        else:
            assert box_length == box_size_i[0], "error, box size is changing"

        bonded_particles, particle_positions, diffusivity = bonds_analysis(trajectory,frame_id=int(frame),step_time=step_time/1e6)

        #skip frames before any bonding
        if(len(bonded_particles)==0):
            continue

        timestep_list.append(step)
        step_time_seconds = step*deltaT_seconds
        time_list.append(step_time_seconds)
        
        #Clustering done below:

        cluster_bondtable,bodies = get_bond_table(trajectory,frame_id=int(frame))
        all_bodies=list(np.arange(0,np.max(bodies)+1))

        clustering = get_clustering_fast(cluster_bondtable)
        combined_list=sorted(flat_list(clustering))

        singletons=find_singletons(combined_list,all_bodies)
        number_singletons=len(singletons)
        singleton_clusters=list(np.ones(number_singletons).astype(int))
            
        clustering=sorted(clustering,key=len)
        cluster_sizes=sorted([len(c) for c in clustering])

        if(min_cluster_size==1):
            cluster_sizes_modified=singleton_clusters+cluster_sizes #add singletons to the list of cluster sizes too only if minimum cluster size chosen is 1
            cluster_sizes=cluster_sizes_modified
        else:
            cluster_sizes_modified=remove_from_list(cluster_sizes,min_cluster_size)  #remove those cluster sizes which are less than the chosen minimum cluster size
            cluster_sizes=cluster_sizes_modified


        if verbose:
            clustering_modified=[[s] for s in singletons]+clustering
            clustering=clustering_modified
            largestclustersizes.append(max(cluster_sizes))
            
        diffusivities_allclusters=[]
        for clus in clustering:
            diffusivities_eachcluster=[]
            for a in clus:
                diff=diffusivity[a]
                diffusivities_eachcluster.append(diff)
            diffusivities_allclusters.append(diffusivities_eachcluster)

        
        mean_diffusivities=[]
        std_diffusivities=[]
        for d in diffusivities_allclusters:
            mean_d=np.mean(d)
            std_d=np.std(d)
            mean_diffusivities.append(mean_d)
            std_diffusivities.append(std_d)

        cluster_sizes_allframes.append(cluster_sizes)
        mean_diffusivities_allframes.append(mean_diffusivities)            
        std_diffusivities_allframes.append(std_diffusivities)

    #Join the data for all frames together into one big list
    
    cluster_sizes_allcombined=[]
    mean_diffusivities_allcombined=[]
    for k in range(len(cluster_sizes_allframes)):
        cluster_sizes=cluster_sizes_allframes[k]
        mean_diffusivities=mean_diffusivities_allframes[k]
        cluster_sizes_allcombined.append(cluster_sizes)
        mean_diffusivities_allcombined.append(mean_diffusivities)
    cluster_sizes_allcombined=flat_list(cluster_sizes_allcombined)
    mean_diffusivities_allcombined=flat_list(mean_diffusivities_allcombined)

    #Grouping by cluster size for the line plot 

    cluster_sizes_allcombined_sorted=sorted(cluster_sizes_allcombined)
    mean_diffusivities_allcombined_sorted = [md for _, md in sorted(zip(cluster_sizes_allcombined, mean_diffusivities_allcombined))]

    grouped_cluster_sizes_allcombined_sorted = group_elements(cluster_sizes_allcombined_sorted)
    grouped_mean_diffusivities_allcombined_sorted = group_elements_basedonalist(cluster_sizes_allcombined_sorted,mean_diffusivities_allcombined_sorted)

    mean_diffusivities_allcombined_sorted_eachclustersize=[]
    std_diffusivities_allcombined_sorted_eachclustersize=[]    
    
    index=0
    cluster_sizes_present=[]

    for entry in grouped_mean_diffusivities_allcombined_sorted:
        mean_value=np.mean(entry)
        std_value=np.std(entry)

        cluster_size_entry=grouped_cluster_sizes_allcombined_sorted[index]
        cluster_sizes_present.append(cluster_size_entry[0])

        mean_diffusivities_allcombined_sorted_eachclustersize.append(mean_value)
        std_diffusivities_allcombined_sorted_eachclustersize.append(std_value)

        index+=1

    return cluster_sizes_allframes,mean_diffusivities_allframes,cluster_sizes_present,mean_diffusivities_allcombined_sorted_eachclustersize,std_diffusivities_allcombined_sorted_eachclustersize

    
if __name__ == "__main__":
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument("--trajectory_file",type=str,required=True)
    parser.add_argument("--min_cluster_size",type=int,required=True)
    parser.add_argument("--minframe",default=1,type=int)
    parser.add_argument("--frameinterval",default=1,type=int)
    parser.add_argument("--cutoff_frames",default=2400,type=int)
    parser.add_argument('--verbose',default=False,action='store_true')

    args = parser.parse_args()
    locals().update(vars(args))
    
    trajectory = gsd.open(trajectory_file,'rb') # read gsd file
    if verbose:
        print('File:',trajectory_file)
        print('Total number of frames:',len(trajectory))

    cluster_sizes_allframes,mean_diffusivities_allframes,cluster_sizes_present,mean_diffusivities_allcombined_sorted_eachclustersize,std_diffusivities_allcombined_sorted_eachclustersize = get_trajectoryinfo_vs_time(trajectory,cutoff_frames,minframe,frameinterval,verbose)

    diffusivity_scatterplotdata=np.array([cluster_sizes_allframes,mean_diffusivities_allframes],dtype=object)
    diffusivity_lineplotdata=np.array([cluster_sizes_present,mean_diffusivities_allcombined_sorted_eachclustersize,std_diffusivities_allcombined_sorted_eachclustersize],dtype=object)

    if(frameinterval==1):
        diffusivity_file1='./../diffusivity_data/skip1frames/'+os.path.splitext(os.path.basename(trajectory_file))[0]+'.diffusivityvsclustersize_minclustersize'+str(min_cluster_size)+'.scatterplot.data'
        diffusivity_scatterplotdata.dump(diffusivity_file1)
        diffusivity_file2='./../diffusivity_data/skip1frames/'+os.path.splitext(os.path.basename(trajectory_file))[0]+'.diffusivityvsclustersize_minclustersize'+str(min_cluster_size)+'.lineplot.data'
        diffusivity_lineplotdata.dump(diffusivity_file2)
    elif(frameinterval==5):
        diffusivity_file1='./../diffusivity_data/skip5frames/'+os.path.splitext(os.path.basename(trajectory_file))[0]+'.diffusivityvsclustersize_minclustersize'+str(min_cluster_size)+'.scatterplot.data'
        diffusivity_scatterplotdata.dump(diffusivity_file1)
        diffusivity_file2='./../diffusivity_data/skip5frames/'+os.path.splitext(os.path.basename(trajectory_file))[0]+'.diffusivityvsclustersize_minclustersize'+str(min_cluster_size)+'.lineplot.data'
        diffusivity_lineplotdata.dump(diffusivity_file2)

    
    
    
    
    
    

