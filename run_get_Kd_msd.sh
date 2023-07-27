#!/bin/bash
nsteps=10000000
#note gamma scale 0.001

submit_binding_simulation () {
    seed=1
    box_length=$1
    volume_fraction_ribosome=$2
    number_rods=$3
    number_linkers=$4
    koff=$5
    binding_distance=$6
    number_gems=$7
    ljeps=0

    otheroptions=""
    final_steps=$nsteps
 
    prefix=prod_v2.6_getKd_newdyn_2023/l${box_length}_vfc${volume_fraction_ribosome}/special_binding_l${box_length}_vfc${volume_fraction_ribosome}_nR${number_rods}_nL${number_linkers}_nG${number_gems}_koff${koff}_lj${ljeps}_bd${binding_distance}_s${seed}

    outprefix=${prefix}_N${final_steps}

    jobname=$(basename $outprefix)
    jobdir=$(dirname $outprefix)
    mkdir -p $jobdir
    simulation_options="--gpu --nsteps $nsteps --koff $koff --box_length $box_length --number_rods $number_rods --number_linkers $number_linkers  --lj_epsilon $ljeps --binding_distance $binding_distance --dump_frequency 20 "

    sbatch --job-name=$jobname -o ${outprefix}.slurm_%j.log -o ${outprefix}.slurm_%j.log --export=outprefix=$outprefix,simulation_options="${simulation_options}" run_simulation_greene_newhoomd_binding.sbatch

}

submit_simulation () {
    seed=1  #make 3 seeds each condition
    box_length=$1
    volume_fraction_ribosome=$2
    number_rods=$3
    number_linkers=$4
    koff=$5
    binding_distance=$6
    number_gems=$7
    ljeps=0

    otheroptions=""
    final_steps=$nsteps
 
    prefix=prod_v2.6_getKd_newdyn_2023/l${box_length}_vfc${volume_fraction_ribosome}/gel_l${box_length}_vfc${volume_fraction_ribosome}_nR${number_rods}_nL${number_linkers}_nG${number_gems}_koff${koff}_lj${ljeps}_bd${binding_distance}_s${seed}

    outprefix=${prefix}_N${final_steps}

    jobname=$(basename $outprefix)
    jobdir=$(dirname $outprefix)
    mkdir -p $jobdir
    simulation_options="--gpu --nsteps $nsteps --koff $koff --box_length $box_length --number_rods $number_rods --number_linkers $number_linkers --volume_fraction_ribosome $volume_fraction_ribosome --lj_epsilon $ljeps --binding_distance $binding_distance --number_gems $number_gems --dump_frequency 20 "

    sbatch --job-name=$jobname -o ${outprefix}.slurm_%j.log -o ${outprefix}.slurm_%j.log --export=outprefix=$outprefix,simulation_options="${simulation_options}" run_simulation_greene_newhoomd.sbatch

}
#for Kd
box_length=400
n_gems=0
for vfc in 0 0.1;do
    for nr in 200;do
         for nl in 200;do
            for koff in 0.001;do
		 for binding_distance in 1.0;do
                    submit_binding_simulation $box_length $vfc $nr $nl $koff $binding_distance $n_gems
		 done
             done
         done
    done
done
exit

#for MSD
box_length=2000
for vfc in 0;do
    for nr in 0;do
         for nl in 0;do
             for koff in 0.0;do
		 for binding_distance in 1.0;do
                     submit_simulation $box_length $vfc $nr $nl $koff $binding_distance 100
		 done
            done
        done
    done
done
exit



#


#langevin tau=3.75e-5 sec
#single step = dt*tau = 2e-3*3.75e-5 = 7.5e-8
#kon = 1/(dt*10) = 1/( 0.002* 10) = 50 /tau = 1.4e6 /sec
#ln(kon/koff) = 
#     ln(50/0.01) = 8.5
#     ln(50/0.001) = 10.8
#     ln(50/0.0001) = 13.1
#     ln(50/0.00001) = 15.4

exit

