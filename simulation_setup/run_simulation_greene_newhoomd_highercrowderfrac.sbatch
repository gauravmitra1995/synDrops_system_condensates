#!/bin/bash
#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --tasks-per-node 1
#SBATCH --gres=gpu:1
#SBATCH --mem 10GB
#SBATCH -t 48:00:00
#SBATCH --dependency=singleton

if [ -z "$outprefix" ];then
    echo "Output prefix required"
    exit
fi
if [ -e "${outprefix}.gsd" ];then
    echo "Skipping $outprefix, already ran"
    exit
fi

set -x
echo "Running in dir: $PWD"
outdir=$(dirname $outprefix)
mkdir -p $outdir

#Regular sim scripts
exe=multivalent_dynbond_v2.65_soft_langevin_compress.py    #to run for VFR 0.4 and 0.5 (any thing above 0.35)

#Kd sim scripts
#exe=multivalent_dynbond_v2.6_getKd_setup_compress.py #to run for VFR 0.4 and 0.5 (any thing above 0.35)

full_outprefix=$(readlink -f $outprefix)

#module purge
#module --ignore-cache load python/intel/3.8.6
#module load hoomd/cuda/openmpi/intel/2.9.3

#python $exe --outprefix $full_outprefix ${simulation_options} 

wrapper=/scratch/projects/hockygroup/data-share/gm2535/pyColloidomer_2023/dybond/run-hoomd2.9.6.bash

if [[ $host =~ ^gm ]]; then
    wrapper=/scratch/projects/hockygroup/data-share/gm2535/pyColloidomer_2023/dybond/run-hoomd2.9.2.bash
fi

$wrapper python $exe --outprefix $full_outprefix ${simulation_options} 

