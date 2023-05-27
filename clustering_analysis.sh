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
current_dir=$(pwd)

wrapper=/scratch/projects/hockygroup/data-share/gm2535/pyColloidomer_2023/dybond/run-hoomd2.9.6.bash

$wrapper python -u get_histograms_clustersize.py --koff $koff --crowder_temperature $crowder_temperature --volume_fraction_ribosome $vfr



 
