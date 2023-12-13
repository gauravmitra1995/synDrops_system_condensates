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

dt=0.002
ng=20
koff0=0
gamma_scale=0.001
sphere_repulsion=500
binding_distance=1.0
cutoff_frames=2400

current_dir=$(pwd)

#wrapper=/scratch/projects/hockygroup/data-share/gm2535/pyColloidomer_2023/dybond/run-hoomd2.9.6.bash


#generate all plots with seeds combined
#python -u combine_data_multipleseeds_and_plot.py --volume_fraction_ribosome $vfr --crowder_temperature $crowder_temperature --koff $koff --min_cluster_size $min_cluster_size --cutoff_frames $cutoff_frames


#largestclustersize vs epsilon
#python -u largestclustersize_vs_epsilon.py --volume_fraction_ribosome $vfr --crowder_temperature $crowder_temperature --koff $koff --min_cluster_size $min_cluster_size --cutoff_frames $cutoff_frames


#largestclustersize vs hexamer conc
python -u largestclustersize_vs_hexamerconc.py --volume_fraction_ribosome $vfr --crowder_temperature $crowder_temperature --number_rods $nr --number_linkers $nl --koff $koff --min_cluster_size $min_cluster_size --cutoff_frames $cutoff_frames

