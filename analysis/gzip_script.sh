#!/bin/bash

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --tasks-per-node 1
#SBATCH --mem 20GB
#SBATCH -t 24:00:00

#parameters to vary:

box_length=$1
shift

vfr=$1
shift

vfp=$1
shift

crowder_temperature=$1
shift

#find /scratch/projects/hockygroup/data-share/gm2535/liam_system/trajectory_data/prod_v2.6_getKd_newdyn_2023/l${box_length}_vfr${vfr}_vfp${vfp}_Tc${crowder_temperature}/*.gsd | wc -l

gzip /scratch/projects/hockygroup/data-share/gm2535/liam_system/trajectory_data/prod_v2.6_getKd_newdyn_2023/l${box_length}_vfr${vfr}_vfp${vfp}_Tc${crowder_temperature}/*.gsd
