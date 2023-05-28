#!/bin/bash
nsteps=100000000
nsteps_original=$nsteps

submit_simulation () {
    koff0=0
    #number_gems=0
    number_gems=20
    gamma_scale=0.001
    box_length=$1
    volume_fraction_ribosome=$2
    volume_fraction_polysome=$3
    number_rods=$4
    number_linkers=$5
    koff=$6
    crowder_temperature=$7
    seed=$8
    binding_distance=1.0
    sphere_repulsion=500
    dt=0.002
    ljeps=0

    if (( $(echo "$volume_fraction_ribosome > 0" |bc -l) ));then
	    walltime="-t 48:00:00"
    fi

    if (( $(echo "$volume_fraction_ribosome > 0.35" |bc -l) ));then
	    nsteps=$(($nsteps_original/2))
	    initial_box_length=1400
	    initial_string=_ls${initial_box_length}
	    initial_command="--initial_box_length $initial_box_length"
    fi

    otheroptions=""
 
    output_prefix=prod_v2.6_newdyn_2023/l${box_length}_vfr${volume_fraction_ribosome}_vfp${volume_fraction_polysome}_Tc${crowder_temperature}/gel_l${box_length}${initial_string}_vfr${volume_fraction_ribosome}_vfp${volume_fraction_polysome}_nG${number_gems}_nR${number_rods}_nL${number_linkers}_k0${koff0}_koff${koff}_repuls${sphere_repulsion}_bd${binding_distance}_Tc${crowder_temperature}_s${seed}_dt${dt}_gs${gamma_scale}

    jobname=$(basename $output_prefix)
    jobdir=$(dirname $output_prefix)
    mkdir -p $jobdir

    restart_simulation_options="--gpu  --nsteps $nsteps --koff $koff --lj_epsilon $ljeps --rod_repulsion $sphere_repulsion --sphere_repulsion $sphere_repulsion --binding_distance $binding_distance --crowder_temperature ${crowder_temperature} --dt $dt --gamma_scale $gamma_scale" #--closed_box"
    simulation_options="--gpu $initial_command --nsteps $nsteps --koff $koff --lj_epsilon $ljeps --binding_distance $binding_distance --rod_repulsion $sphere_repulsion --sphere_repulsion $sphere_repulsion --crowder_temperature ${crowder_temperature} --number_rods $number_rods --number_linkers $number_linkers --volume_fraction_ribosome $vfr --volume_fraction_polysome $vfp --number_gems $number_gems --dt=$dt --gamma_scale $gamma_scale" #--closed_box"

    sbatch $walltime --job-name=$jobname -o ${output_prefix}.slurm_%j.log --export=outprefix=$output_prefix,nsteps=$nsteps,simulation_options="${simulation_options}",restart_simulation_options="${restart_simulation_options}"  run_simulation_greene_newhoomd_restart.sbatch

}

box_length=860
#box_length=400
#nr = 1170, nl=390 are default
nr=1170
nl=390
#for vfr_vfp in 0.0_0 0.15_0 0.2_0 0.25_0 0.3_0 0.35_0 0.4_0 0.5_0;do
for vfr_vfp in 0.3_0;do
#for vfr_vfp in 0.4_0 0.5_0;do #0.5_0;do
#for vfr_vfp in 0.4_0 0.5_0;do
#for vfr_vfp in 0.15_0 0.25_0;do
    #for koff in 0.001 0.0003 0.0001;do
    for koff in 0.001;do
    #for koff in 0.0001;do
        vfr=$(echo $vfr_vfp |cut -f 1 -d '_')
        vfp=$(echo $vfr_vfp |cut -f 2 -d '_')
	#for crowder_temperature in 1.0;do
	#for crowder_temperature in 0.5 2.0;do
	for crowder_temperature in 1.1 1.2;do
            for seed in `seq 1 5`;do 
                 submit_simulation $box_length $vfr $vfp $nr $nl $koff $crowder_temperature $seed
            done
        done
    done
done
