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

"""
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


"""
time_convert=7.5e-8
for frame in frames_list:
    system=trajectory[int(frame)]
"""


