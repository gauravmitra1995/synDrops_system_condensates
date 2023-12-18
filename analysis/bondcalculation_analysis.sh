#!/bin/bash

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 2GB
#SBATCH -t 1:00:00
##SBATCH --dependency=singleton

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


dt=0.002
ng=0
koff0=0
gamma_scale=0.001
sphere_repulsion=500
binding_distance=1.0
 
kon=50.0

if [ $kon -eq 0 ];then
       epsilon=0.0
else
      if [ $koff -eq 0 ];then
             epsilon='infinite'
      else
             epsilon_raw=$(echo "l($kon/$koff)" | bc -l)
             epsilon=`printf "%.1f" $epsilon_raw`
      fi
fi

current_dir=$(pwd)

#wrapper=/scratch/projects/hockygroup/data-share/gm2535/pyColloidomer_2023/dybond/run-hoomd2.9.6.bash

for seed in `seq 1 5`;do 

      python -u bondcount_Kd_calculation.py --trajectory_file ./gsdfiles_getKd/gel_l${box_length}_vfr${vfr}_vfp${vfp}_nG${ng}_nR${nr}_nL${nl}_k0${koff0}_koff${koff}_repuls${sphere_repulsion}_bd${binding_distance}_Tc${crowder_temperature}_s${seed}_dt${dt}_gs${gamma_scale}.allruns.gsd --epsilon $epsilon

done
