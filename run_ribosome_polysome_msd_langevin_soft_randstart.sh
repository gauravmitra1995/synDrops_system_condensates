#!/bin/bash
nsteps=100000000

submit_simulation () {
    koff0=0
    #number_gems=0
    number_gems=20
    gamma_scale=0.001
    prev_steps=0
    #prev_steps=$nsteps
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

    otheroptions=""
    final_steps=$(($prev_steps+$nsteps))
 
    input_prefix=prod_v2.6_newdyn_2023/l${box_length}_vfr${volume_fraction_ribosome}_vfp${volume_fraction_polysome}_Tc${crowder_temperature}/gel_l${box_length}_vfr${volume_fraction_ribosome}_vfp${volume_fraction_polysome}_nG${number_gems}_nR${number_rods}_nL${number_linkers}_k0${koff0}_koff${koff0}_repuls${sphere_repulsion}_bd${binding_distance}_Tc${crowder_temperature}_s${seed}_dt${dt}_gs${gamma_scale}_N${prev_steps}
    output_prefix=prod_v2.6_newdyn_2023/l${box_length}_vfr${volume_fraction_ribosome}_vfp${volume_fraction_polysome}_Tc${crowder_temperature}/gel_l${box_length}_vfr${volume_fraction_ribosome}_vfp${volume_fraction_polysome}_nG${number_gems}_nR${number_rods}_nL${number_linkers}_k0${koff0}_koff${koff}_repuls${sphere_repulsion}_bd${binding_distance}_Tc${crowder_temperature}_s${seed}_dt${dt}_gs${gamma_scale}

    outprefix=${output_prefix}_N${final_steps}

    jobname=$(basename $outprefix)
    jobdir=$(dirname $outprefix)
    mkdir -p $jobdir

    if [ -e "${input_prefix}.gsd" ];then
        simulation_options="--gpu -i ${input_prefix}.gsd --nsteps $nsteps --koff $koff --lj_epsilon $ljeps --rod_repulsion $sphere_repulsion --sphere_repulsion $sphere_repulsion --binding_distance $binding_distance --crowder_temperature ${crowder_temperature} --dt $dt --gamma_scale $gamma_scale" #--closed_box"

        sbatch $walltime --job-name=$jobname -o ${outprefix}.slurm_%j.log -o ${outprefix}.slurm_%j.log --export=outprefix=$outprefix,simulation_options="${simulation_options}"  run_simulation_greene_newhoomd.sbatch
	#run_simulation_greene.sbatch
    else
        simulation_options="--gpu --nsteps $nsteps --koff $koff --lj_epsilon $ljeps --binding_distance $binding_distance --rod_repulsion $sphere_repulsion --sphere_repulsion $sphere_repulsion --crowder_temperature ${crowder_temperature} --number_rods $number_rods --number_linkers $number_linkers --volume_fraction_ribosome $vfr --volume_fraction_polysome $vfp --number_gems $number_gems --dt=$dt --gamma_scale $gamma_scale" #--closed_box"
        sbatch $walltime --job-name=$jobname -o ${outprefix}.slurm_%j.log -o ${outprefix}.slurm_%j.log --export=outprefix=$outprefix,simulation_options="${simulation_options}"  run_simulation_greene_newhoomd.sbatch
    fi


}

box_length=860
#box_length=400
#nr = 1170, nl=390 are default
nr=1170
nl=390
#for vfr_vfp in 0.0_0 0.15_0 0.2_0 0.25_0 0.3_0 0.35_0;do
for vfr_vfp in 0.3_0;do
    #for koff in 0.001 0.0003 0.0001;do
    for koff in 0.001;do
        vfr=$(echo $vfr_vfp |cut -f 1 -d '_')
        vfp=$(echo $vfr_vfp |cut -f 2 -d '_')
	#for crowder_temperature in 1.0;do
	for crowder_temperature in 0.5 0.9 1.0 1.1 1.5 2.0;do
            for seed in `seq 1 3`;do 
                 submit_simulation $box_length $vfr $vfp $nr $nl $koff $crowder_temperature $seed
            done
        done
    done
done

exit

#vary koff
box_length=860
#box_length=400
#nr = 1170, nl=390 are default
nr=1170
nl=390
crowder_temperature=0.9
#koff=0
for vfr_vfp in 0.15_0 0.2_0 0.3_0;do
#for vfr_vfp in 0.0_0;do
    vfr=$(echo $vfr_vfp |cut -f 1 -d '_')
    vfp=$(echo $vfr_vfp |cut -f 2 -d '_')
    #for koff in 0.0000 0.0001 0.0002 0.0003 0.0004 0.0005 0.0008;do
    for koff in 0.0003;do
            for binding_distance in 0.9;do
                 submit_simulation $box_length $vfr $vfp $nr $nl $koff $crowder_temperature $binding_distance
            done
    done
done
exit

exit


#small
box_length=596.3
#box_length=400
#nr = 1170, nl=390 are default
nr=390
nl=130
koff=0.0002
#koff=0
for vfr_vfp in 0.0_0 0.15_0 0.3_0;do
    vfr=$(echo $vfr_vfp |cut -f 1 -d '_')
    vfp=$(echo $vfr_vfp |cut -f 2 -d '_')
	 for crowder_temperature in 0.9 1.0 1.1;do
            for binding_distance in 0.9;do
                 submit_simulation $box_length $vfr $vfp $nr $nl $koff $crowder_temperature $binding_distance
            done
        done
done
exit


exit


exit

