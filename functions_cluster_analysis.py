import os
import sys
#import hoomd
import pickle
import numpy as np
from gsd import hoomd as gsd
import argparse
import matplotlib
#matplotlib.use('svg')
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import json
import pandas as pd
from numpy.lib.stride_tricks import sliding_window_view
sns.set_context("poster", font_scale=1.0)

def get_clustering(cluster_bondtable_new):
    topo=cluster_bondtable_new
    l=len(topo)
    q=[]
    q=topo

    output = []
    while len(q)>0:
        first, *rest = q
        first = set(first)

        lf = -1
        while len(first)>lf:
            lf = len(first)

            rest2 = []
            for r in rest:
                if len(first.intersection(set(r)))>0:
                    first |= set(r)
                else:
                    rest2.append(r)
            rest = rest2

        output.append(sorted(list(first)))
        q = rest

    return(output)

def window_average(values, window_size, stride=None, xvalues=None):
    """
        This function returns a window average of data over length `window_size', and returns the value every `stride' steps. 
        If xvalues is given, it also returns the window averaged x-values, which should be the location of the averaged points.
    """
    if stride is None:
        stride=window_size
    window_avg = sliding_window_view(values,window_size).mean(axis=1)[::stride]
    window_std = sliding_window_view(values,window_size).std(axis=1)[::stride]
    if xvalues is not None:
        window_avg_x = sliding_window_view(xvalues,window_size).mean(axis=1)[::stride]
    else:
        window_avg_x = None
    return window_avg, window_std, window_avg_x


def find_cluster1(Alist,tag): # given tag, find which cluster it belongs to, outputs cluster number
    if tag>=Alist[-1]:
        return(len(Alist)-1)
    else:
        for n in range(len(Alist)):
            if tag < Alist[n]:
                return(n-1)
def remove_repeats(listt): # [[1,2],[1,2],[2,1],[4,3]] --> [[1,2],[4,3]]
    result=list()
    countrepeats={}
    for i in listt:
        countrepeats[str(min(i[0],i[1]))+','+str(max(i[0],i[1]))]=0
        if i not in result and i[::-1] not in result:
            result.append(i)
        elif i in result or i[::-1] in result:
            countrepeats[str(min(i[0],i[1]))+','+str(max(i[0],i[1]))]+=1
    return(result,list(countrepeats.values()))
    #return result

def count_repeats_particles(listt): 
    for i in range(len(listt)):
        particlei=listt[i][0]
        for j in range(i,len(listt)):
            particlej0=listt[j][0]
            particlej1=listt[j][1]
            if(i!=j):
                if(particlei==particlej0 or particlei==particlej1):
                    print('Found a match:')
                    print(listt[i],listt[j]) 
                    print("****************************************************")       


def flat_list(alist): # [[1,2],[3,4]] --> [1,2,3,4]
    flatlist=list()
    for sublist in alist:
        for i in sublist:
            flatlist.append(i)
    return(flatlist)

def makepairs(listt): #[1,2,3]-->[[1,2],[1,3],[2,3]]
    result=list()
    for i in range(len(listt)):
        for j in range(len(listt)):
            if i < j:
                result.append([listt[i],listt[j]])
    return(result)

def makepairs2(listt):  
    if len(listt)==2: # [1,2]-->[1,2]
        return([listt])

    else: #[1,2,3]-->[[1,2],[2,3],[3,1]];
        result=list()
        for i in range(len(listt)):
            result.append([listt[i],listt[(i+1)%len(listt)]])
        return(result)

def counting_dict(listt):  # [1,1,0,2,2,2] --> {1:2, 0:1, 2:3}
    result={}
    for i in listt:
        if i in result.keys():
            result[i]+=1
        else:
            result[i]=1
    return(result)
