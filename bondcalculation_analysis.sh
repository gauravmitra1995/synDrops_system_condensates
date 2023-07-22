#!/bin/bash

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

current_dir=$(pwd)

#wrapper=/scratch/projects/hockygroup/data-share/gm2535/pyColloidomer_2023/dybond/run-hoomd2.9.6.bash

for seed in `seq 1 5`;do 

      python -u bondcount_Kd_calculation.py --trajectory_file ./gsdfiles_getKd/gel_l${box_length}_vfr${vfr}_vfp${vfp}_nG${ng}_nR${nr}_nL${nl}_k0${koff0}_koff${koff}_repuls${sphere_repulsion}_bd${binding_distance}_Tc${crowder_temperature}_s${seed}_dt${dt}_gs${gamma_scale}.allruns.gsd

done
