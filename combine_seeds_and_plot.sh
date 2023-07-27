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

min_cluster_size=$1
shift

#cutoff_frames=2400

dt=0.002
ng=20
koff0=0
gamma_scale=0.001
sphere_repulsion=500
binding_distance=1.0

current_dir=$(pwd)

#wrapper=/scratch/projects/hockygroup/data-share/gm2535/pyColloidomer_2023/dybond/run-hoomd2.9.6.bash

python -u combine_data_multipleseeds_and_plot.py --volume_fraction_ribosome $vfr --crowder_temperature $crowder_temperature --koff $koff --min_cluster_size $min_cluster_size 

