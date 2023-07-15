from functions_cluster_analysis import *

parser=argparse.ArgumentParser()
parser.add_argument("--trajectory_file",type=str)
parser.add_argument("--minframe",default=0,type=int)
#parser.add_argument("--maxframe",default=100,type=int)
parser.add_argument("--dt",type=float)
parser.add_argument("--frameinterval",default=1,type=int)
args = parser.parse_args()
locals().update(vars(args))

trajectory = gsd.open(trajectory_file,'rb') # read gsd file
print('file:',trajectory_file)
print('\n')

num_frames=len(trajectory)
print('Total number of frames')
print(num_frames)

maxframe=num_frames-1
"""
######################################## Finding last uncorrupted/readable frame to use as maxframe #############################################

for i in reversed(range(num_frames)):
    print('trying frame',i)
    try:
        last_frame_testpos=trajectory[i].particles.position[0]
        last_frame=i
        break
    except:
        print('next frame')

if last_frame < maxframe:
    print('the given maxframe is unreadable, using the last readable frame, frame',last_frame,'as maxframe')
    maxframe=last_frame
"""

##################################################################################################################################################
frames_list=np.arange(minframe,maxframe+1,frameinterval) # list of frames to analyze

####################################################### Gathering basic info #######################################################

system = trajectory[0] # gather basic information from 0th frame

box_size = list(system.configuration.box[:3])

type_list = list(system.particles.types) # list of all particle types in the system
num_particles = system.particles.N # total num of particles in the system

system0 = trajectory[0]
print("Timestep at frame 0:", system0.configuration.step)
system1 = trajectory[1]
print("Timestep at frame 1:", system1.configuration.step)
system2 = trajectory[2]
print("Timestep at frame 2:", system2.configuration.step)
systemfinalminus1=trajectory[-2]
print("Timestep at frame before final:",systemfinalminus1.configuration.step)
systemfinal = trajectory[-1]
print("Timestep at final frame:",systemfinal.configuration.step)


"""
time_convert=7.5e-8
for frame in frames_list:
    system=trajectory[int(frame)]
    if(frame==frames_list[1600]):
        print("Frame no, simulation timestep:")
        print(frame-2,trajectory[int(frame-2)].configuration.step)
        print(frame-1,trajectory[int(frame-1)].configuration.step)
        print(frame,trajectory[int(frame)].configuration.step)
        print(frame+1,trajectory[int(frame+1)].configuration.step)
        print("Time in seconds at frame 1599:")
        print(trajectory[int(frame)-1].configuration.step*time_convert)
        #print("Difference in timesteps between consecutive frames:")
        #print(trajectory[int(frame)].configuration.step-trajectory[int(frame)-1].configuration.step)
        #print("Difference in timesteps between consecutive frames (in units of sec):")
        #print((trajectory[int(frame)].configuration.step-trajectory[int(frame)-1].configuration.step)*time_convert)
"""


"""
bonds_allframes=[]

A_tags_list=[i for i in range(num_particles) if system.particles.typeid[i]==0] # tags of central particles
num_clusters=len(A_tags_list) # number of clusters
Np=int((num_particles/num_clusters-1)/2)
cluster_ids = np.arange(num_clusters,dtype=int)

typeid_list = [list(system.particles.typeid[i : j]) for i, j in zip([None]+A_tags_list, A_tags_list+[None]) if i!=None] # [[..cluster0 typeids..],[..cluster1 typeids..],..]


Np_list=[(len(i)-1)/2 for i in typeid_list] # list of Np's for each cluster
terminal_typeids=[i[-1] for i in typeid_list] # type ids of terminal particles (based on last typeid in each typeid sublist in typeid_list)

terminal_types=[type_list[i] for i in terminal_typeids]
terminal_types_unindexed=[i[0] for i in terminal_types]

total_clusterharmonicbonds=int(sum([i*2 for i in Np_list]))

num_non_dybonds=total_clusterharmonicbonds

bond_types=list(system.bonds.types)

dybond_types=[i for i in system.bonds.types if ((i[0] in terminal_types_unindexed) and (i[1] in terminal_types_unindexed)) or ('-' in i) or (',') in i] # if the bond type name has '-' in it, it's a dybond (maybe use some other criteria)

dybond_typeids=[system.bonds.types.index(i) for i in dybond_types]

"""

