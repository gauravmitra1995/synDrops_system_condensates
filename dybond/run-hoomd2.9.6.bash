#!/bin/bash

scriptdir=$(cd $(dirname $0);pwd)

args=''
for i in "$@"; do
    i="${i//\\/\\\\}"
    args="$args \"${i//\"/\\\"}\""
done

if [ "$args" == "" ]; then args="/bin/bash"; fi

if [ -e /dev/nvidia0 ]; then nv="--nv"; fi

singularity \
    exec $nv \
    --bind $scriptdir/nvidia-hoomd-2.9.6/update.py:/ext3/hoomd/hoomd/dybond_plugin/update.py:ro \
    --bind $scriptdir/nvidia-hoomd-2.9.6/_dybond_plugin.so:/ext3/hoomd/hoomd/dybond_plugin/_dybond_plugin.so:ro \
    --overlay $scriptdir/singularity/hoomd-2.9.6.sqf:ro \
    --overlay $scriptdir/singularity/hoomd2.9.6-dep-20210524.sqf:ro \
    --overlay $scriptdir/singularity/hoomd-2.9.6-20220414.ext3:ro \
    $scriptdir/singularity/cuda11.1.1-cudnn8-devel-ubuntu20.04.sif \
    /bin/bash -c "
source /ext3/env.sh
$args
exit
"
