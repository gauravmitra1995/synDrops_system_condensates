#!/bin/bash
#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --tasks-per-node 1
#SBATCH --gres=gpu:1
#SBATCH --mem 20GB
#SBATCH -t 24:00:00
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

exe=get_binding_constant_v2.3.py
full_outprefix=$(readlink -f $outprefix)

#module purge
#module --ignore-cache load python/intel/3.8.6
#module load hoomd/cuda/openmpi/intel/2.9.3

#export PYTHONPATH=/scratch/work/hockygroup/software/hoomd-v2.9.7:$PYTHONPATH
#export PYTHONPATH=/scratch/work/hockygroup/software/hoomd-v2.9.7-dynbond/dybond_plugin:$PYTHONPATH
#python $exe --outprefix $full_outprefix ${simulation_options} 

/scratch/projects/hockygroup/data-share/gm2535/pyColloidomer_2023/dybond/run-hoomd2.9.6.bash python $exe --outprefix $full_outprefix ${simulation_options} 

