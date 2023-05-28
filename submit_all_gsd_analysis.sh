#!/bin/bash
#SBATCH --mem 20GB
#SBATCH --time 08:00:00

module load python/intel/3.8.6
source ~/pyenv/md-py3.8/bin/activate

out_dir=prod_v2.6_newdyn_2023
for dir in $(find ${out_dir} -maxdepth 1 -mindepth 1 -type d);do 
	echo $dir
	gsdfiles=$(ls $dir/*.gsd| grep koff0.001 |grep -e Tc0.5 -e Tc1.0 -e Tc2.0)
	prefixes=$(echo $gsdfiles |xargs -n 1 basename |cut -f 1-15 -d '_' |sort | uniq )
	for prefix in $prefixes;do
		gsdfiles=$(ls -rt $dir/${prefix}*.gsd |grep -v compress |grep -v combined)
		nfiles=$(echo $gsdfiles |wc -w)
		if [ $nfiles -gt 1 ];then
			newfile=$dir/${prefix}_combined.gsd
		        python ./concatenate_gsd.py --stride 5 --force $newfile $gsdfiles
		        bash submit_gsd_analysis.sh $newfile
		else
		    bash submit_gsd_analysis.sh $gsdfiles
		fi
	done
done

