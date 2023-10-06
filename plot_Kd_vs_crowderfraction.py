from functions_cluster_analysis import *

def calculate_Kd(bonds,box_length,Navo,nr):
    volume=(box_length*(10**-8))**3
    Kd= ((nr-bonds)/(Navo*volume))*(float(nr/bonds)-1.0)
    return Kd

parser=argparse.ArgumentParser()
#parser.add_argument("--vfr",default=0.0,type=float)
parser.add_argument("--epsilon",default=10.8,required=True,type=float)
parser.add_argument("--koff",default=0.001,required=True,type=float)
parser.add_argument("--box_length",default=400,type=int)
parser.add_argument("--nr",default=200,type=int)
parser.add_argument("--nl",default=200,type=int)
#parser.add_argument("--seed",default=1,type=int)

args = parser.parse_args()
locals().update(vars(args))



# Define the directory path and the pattern
directory_path = './Kd_data/'
pattern = 'gel_l'+str(box_length)+'_vfr*_vfp0_nG0_nR'+str(nr)+'_nL'+str(nl)+'_k00_koff'+str(koff)+'_repuls500_bd1.0_Tc1.0_s*_dt0.002_gs0.001.allruns.bondsatsaturation.data'

# Use glob to find all files matching the pattern
matching_files = glob.glob(directory_path + pattern)

# Define a custom key function to extract 's' and 'vfr' values for sorting
def custom_sort_key(file_path):
    # Extract the value of 's' from the file path (assuming it's the last integer before ".allruns.bondsatsaturation.data")
    s = int(file_path.split("_s")[-1].split("_")[0])

    # Extract the value of 'vfr' from the file path (assuming it's the first floating-point number after "vfr")
    vfr = float(file_path.split("vfr")[1].split("_")[0])

    return (vfr, s)

# Sort the file paths using the custom key function
sorted_files = sorted(matching_files, key=custom_sort_key)

# Create a dictionary to store the combined data for each 'vfr' value
combined_data = {}

# Read the data from each file and combine it based on 'vfr' values
for file_path in sorted_files:
    data = np.load(file_path, allow_pickle=True)

    # Extract 'vfr' from the file path
    vfr = float(file_path.split("vfr")[1].split("_")[0])
    #print(vfr,data[2])

    # Append the data to the corresponding 'vfr' entry in the dictionary
    if vfr not in combined_data:
        combined_data[vfr] = [data[2]]
    else:
        combined_data[vfr].append(data[2])

Navo=6.023e23
# Print the combined data for each 'vfr' value
vfr_values=[]
mean_Kd_values=[]
stddev_Kd_values=[]

Kd_values_micromolar_allcrowderfractions={}

for vfr, bonds in combined_data.items():
    print(f"Bond count for vfr={vfr}:")
    bonds=np.array(bonds)
    print("Bonds: ",bonds)
    Kd_values=np.zeros(bonds.shape[0])
    for i in range(bonds.shape[0]):
        b=bonds[i]
        Kd_values[i]=calculate_Kd(b,box_length,Navo,nr)

    Kd_values_micromolar=Kd_values*10**6
    Kd_values_micromolar_allcrowderfractions[vfr]=Kd_values_micromolar

    print("Kd values (in micro molar): ",Kd_values_micromolar)
    print("-"*100)

    vfr_values.append(vfr)
    mean_Kd_values.append(np.mean(Kd_values_micromolar,axis=0))
    stddev_Kd_values.append(np.std(Kd_values_micromolar,axis=0))


#print(Kd_values_micromolar_allcrowderfractions)

"""
# Print the sorted file paths
for file_path in sorted_files:
    s = int(file_path.split("_s")[-1].split("_")[0])
    vfr = float(file_path.split("vfr")[1].split("_")[0])
    #print(vfr,s)
"""
#colors=['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd']

fig=plt.figure(figsize=(8,6),dpi=200)
plt.errorbar(vfr_values,mean_Kd_values,yerr=stddev_Kd_values,xerr=None,marker='o',markersize=6.0,linewidth=1.0,color='black',elinewidth=1.0,capsize=5.0,ecolor='red',capthick=1.0)
r=0
for vfr in vfr_values:
    plt.scatter([vfr]*5,Kd_values_micromolar_allcrowderfractions[vfr],marker='^',s=20.0,color='blue')
    r+=1

plt.xlabel(r'$\phi_{crowder}$')
plt.ylabel(r'$K_{D}(\mu M)$')
plt.grid(alpha=0.5)
#plt.yticks(np.arange(0,11,2))
plt.title(r'$\varepsilon = $'+str(epsilon))
fig.tight_layout()
plt.savefig('./final_figures/Kd_vs_crowderfrac/plot_epsilon'+str(epsilon)+'_Kdvscrowderfraction.svg')
plt.close()

