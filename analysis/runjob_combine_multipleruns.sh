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

binding_distance=$1
shift

dt=0.002
#ng=0  #For KD simulations, no of gems was 0
ng=20
koff0=0
gamma_scale=0.001
sphere_repulsion=500

current_dir=$(pwd)

wrapper=/scratch/projects/hockygroup/data-share/gm2535/pyColloidomer_2023/dybond/run-hoomd2.9.6.bash

for seed in `seq 1 5`;do 


     #fileprefix=/scratch/projects/hockygroup/data-share/gm2535/liam_system/prod_v2.6_getKd_newdyn_2023/l${box_length}_vfr${vfr}_vfp${vfp}_Tc${crowder_temperature}/gel_l${box_length}_vfr${vfr}_vfp${vfp}_nG${ng}_nR${nr}_nL${nl}_k0${koff0}_koff${koff}_repuls${sphere_repulsion}_bd${binding_distance}_Tc${crowder_temperature}_s${seed}_dt${dt}_gs${gamma_scale}
     fileprefix=/scratch/work/hockygroup/data-share/gm2535/liam_system/prod_v2.6_newdyn_2023/l${box_length}_vfr${vfr}_vfp${vfp}_Tc${crowder_temperature}/gel_l${box_length}_vfr${vfr}_vfp${vfp}_nG${ng}_nR${nr}_nL${nl}_k0${koff0}_koff${koff}_repuls${sphere_repulsion}_bd${binding_distance}_Tc${crowder_temperature}_s${seed}_dt${dt}_gs${gamma_scale}
      
     if (( $(echo "$vfr > 0.35" |bc -l) ));then 

      #fileprefix=/scratch/projects/hockygroup/data-share/gm2535/liam_system/prod_v2.6_getKd_newdyn_2023/l${box_length}_vfr${vfr}_vfp${vfp}_Tc${crowder_temperature}/gel_l${box_length}_ls1400_vfr${vfr}_vfp${vfp}_nG${ng}_nR${nr}_nL${nl}_k0${koff0}_koff${koff}_repuls${sphere_repulsion}_bd${binding_distance}_Tc${crowder_temperature}_s${seed}_dt${dt}_gs${gamma_scale}
      fileprefix=/scratch/work/hockygroup/data-share/gm2535/liam_system/prod_v2.6_newdyn_2023/l${box_length}_vfr${vfr}_vfp${vfp}_Tc${crowder_temperature}/gel_l${box_length}_ls1400_vfr${vfr}_vfp${vfp}_nG${ng}_nR${nr}_nL${nl}_k0${koff0}_koff${koff}_repuls${sphere_repulsion}_bd${binding_distance}_Tc${crowder_temperature}_s${seed}_dt${dt}_gs${gamma_scale}

     fi

     $wrapper python -u join_gsdfiles_severalruns.py --fileprefix $fileprefix --vfr $vfr

done
