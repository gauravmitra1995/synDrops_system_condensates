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

exe=multivalent_dynbond_v2.65_soft_langevin_compress.py

prefix=$(readlink -f $outprefix)
progressfile=${prefix}.progress.txt
prgtmp=$(cat $progressfile)
if [ -z "$prgtmp" ];then
    prev_steps=0
    final_steps=$nsteps
    run_options=${simulation_options}
else
    prev_steps=$prgtmp
    final_steps=$(($prev_steps+$nsteps))
    input_prefix=${prefix}_N$prev_steps
    extra_flag="-i ${input_prefix}.gsd"
    run_options=${restart_simulation_options}
fi

full_outprefix=${prefix}_N$final_steps

/scratch/projects/hockygroup/data-share/gm2535/pyColloidomer_2023/dybond/run-hoomd2.9.6.bash python $exe $extra_flag --outprefix $full_outprefix ${run_options} $extra_flag && echo $final_steps > $progressfile

