from functions_cluster_analysis import *
from scipy.optimize import curve_fit

# Define the equation B(t) = b (1 - exp(-kt))
def equation(t, b, k):
    return b * (1 - np.exp(-k * t))

parser=argparse.ArgumentParser()
parser.add_argument("--trajectory_file",type=str,required=True)
parser.add_argument("--minframe",default=0,type=int)
#parser.add_argument("--maxframe",default=100,type=int)
#parser.add_argument("--dt",default=0.002,type=float)

parser.add_argument("--frameinterval",default=1,type=int)
#parser.add_argument("--vfr",default=0.0,type=float,required=True)
parser.add_argument("--epsilon",default=10.8,type=float,required=True)
#parser.add_argument("--koff",default=0.001,type=float,required=True)
#parser.add_argument("--seed",default=1,type=int)

args = parser.parse_args()
locals().update(vars(args))

trajectory = gsd.open(trajectory_file,'rb') # read gsd file
print('file:',trajectory_file)
print('\n')
print('Epsilon:',epsilon)

vfr=float(os.path.basename(trajectory_file).split("_")[2].replace('vfr',''))
print("Crowder volume fraction: ",vfr)

seed=int(os.path.basename(trajectory_file).split("_")[12].replace('s',''))
print("Seed: ",seed)

dt=float(os.path.basename(trajectory_file).split("_")[13].replace('dt',''))
print("dt: ",dt)

koff=float(os.path.basename(trajectory_file).split("_")[8].replace('koff',''))
print("koff: ",koff)

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
print("Box size:")
print(box_size)

type_list = list(system.particles.types) # list of all particle types in the system
num_particles = system.particles.N # total num of particles in the system
bonds_allframes=[]
time_allframes=[]
steptime_in_sec=7.5e-2/1e6
timesteps_allframes=[]

for frame in frames_list:
    system=trajectory[int(frame)]
    bonds_allframes.append(system.bonds.N)
    timestep=system.configuration.step
    time_in_sec=steptime_in_sec*timestep
    time_allframes.append(time_in_sec)
    timesteps_allframes.append(timestep)


time_allframes=np.array([0.0]+time_allframes)
timesteps_allframes=np.array([0]+timesteps_allframes).astype(int)
bonds_allframes=np.array([0]+bonds_allframes).astype(int)

params, _ = curve_fit(equation, time_allframes, bonds_allframes)
# Extract the fitted parameters
b_fit, k_fit = params

print("VFR, SEED, FITTED BOND COUNT")
plot_data=np.array([vfr,seed,b_fit])
print(plot_data)

bondcount_data_file='./Kd_data/'+os.path.splitext(os.path.basename(trajectory_file))[0]+'.bondsatsaturation.data'
plot_data.dump(bondcount_data_file)

bonds_fitted=equation(time_allframes, b_fit, k_fit) 

fig=plt.figure(figsize=(8,6),dpi=100)
plt.plot(time_allframes[::2],bonds_allframes[::2],marker='o',markersize=2.5,linewidth=1.0,color='k',label='Data')
plt.plot(time_allframes,bonds_fitted,marker=None,linewidth=2.5,linestyle='--',color='red',label='Fit')
plt.xlabel('Time (in seconds)')
plt.ylabel('# Bonds')
plt.grid(alpha=0.5)
plt.yticks(np.arange(0,210,20))
plt.legend(loc='lower right')
plt.title(r'$\phi_{ribosome} = $'+str(vfr)+' ; '+r'$\varepsilon = $'+str(epsilon)+' ; '+r'$s = $'+str(seed))
fig.tight_layout()
plt.xticks(np.arange(0,1.55,0.3))
plt.savefig('./final_figures/Kd_vs_crowderfrac/plot_vfr'+str(vfr)+'_epsilon'+str(epsilon)+'_seed'+str(seed)+'_bondcountvstime_withfit.svg',bbox_inches='tight')
plt.close()
