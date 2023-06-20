#!/bin/bash

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 4GB
#SBATCH -t 5:00:00
##SBATCH --dependency=singleton
##SBATCH --gres=gpu:1
##SBATCH --gres=gpu:a100:1

#parameters to vary:

box_length=$1
shift

vfr=$1
shift

vfp=$1
shift

nr=$1
shift

nl=$1
shift

koff=$1
shift

crowder_temperature=$1
shift

min_cluster_size=$1
shift

dt=0.002
ng=20
koff0=0
gamma_scale=0.001
sphere_repulsion=500
binding_distance=1.0

current_dir=$(pwd)

#wrapper=/scratch/projects/hockygroup/data-share/gm2535/pyColloidomer_2023/dybond/run-hoomd2.9.6.bash

for seed in `seq 1 5`;do 

      trajectory_file=/scratch/work/hockygroup/gm2535/prod_v2.6_newdyn_2023/l${box_length}_vfr${vfr}_vfp${vfp}_Tc${crowder_temperature}/gel_l${box_length}_vfr${vfr}_vfp${vfp}_nG${ng}_nR${nr}_nL${nl}_k0${koff0}_koff${koff}_repuls${sphere_repulsion}_bd${binding_distance}_Tc${crowder_temperature}_s${seed}_dt${dt}_gs${gamma_scale}_combined.gsd


      if (( $(echo "$vfr > 0.35" |bc -l) ));then

         trajectory_file=/scratch/work/hockygroup/gm2535/prod_v2.6_newdyn_2023/l${box_length}_vfr${vfr}_vfp${vfp}_Tc${crowder_temperature}/gel_l${box_length}_ls1400_vfr${vfr}_vfp${vfp}_nG${ng}_nR${nr}_nL${nl}_k0${koff0}_koff${koff}_repuls${sphere_repulsion}_bd${binding_distance}_Tc${crowder_temperature}_s${seed}_dt${dt}_combined.gsd

      fi


      python -u clustering_using_bondtable_2023_v2.py --trajectory_file $trajectory_file --verbose --min_cluster_size $min_cluster_size

done
