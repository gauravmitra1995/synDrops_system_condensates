#!/bin/bash
gsd_file=$1
if [ ! -e "$gsd_file" ];then
    echo "Usage: $0 gsdfile"
    exit 0
fi

job_name=results_$(basename $gsd_file .gsd)
log_file=$(dirname $gsd_file)/${job_name}.log

sbatch --job-name $job_name -o $log_file --export gsd_file=$gsd_file analyze_gsd_papermill.sh
