#!/bin/bash
#SBATCH --mem=20GB
#SBATCH --nodes=1
#SBATCH --tasks-per-node 1
#SBATCH --cpus-per-task 4
#SBATCH --time 24:00:00
module load python/intel/3.8.6
source ~/pyenv/md-py3.8/bin/activate
analysis_notebook=/scratch/projects/hockygroup/gmh4/projects/hoomd-liam-gel/ab_sim_graphs_for_Glen_v4.51_directgsd_papermill.ipynb

#gsd_file=$1
if [ ! -e "$gsd_file" ];then
    echo "Usage: $0 gsdfile"
    exit 0
fi

run_dir=$(dirname $gsd_file)
output_notebook=$run_dir/results_v4.51_$(basename $gsd_file .gsd).ipynb
#check if output notebook exists if don't want to repeat work
#if [ ! -e "$output_notebook" ];then
papermill $analysis_notebook $output_notebook -r gsd_file $gsd_file
#fi

