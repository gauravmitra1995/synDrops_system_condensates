#!/bin/bash
prefix=$1

max_number=$(ls ${prefix}_N*.gsd|grep -v compress|awk -F'_' '{print $NF}' | tr -d 'N' |tr -d '.gsd'|sort -n |tail -n 1)
if [ ! -z ${max_number} ];then
    echo $max_number > ${prefix}.progress.txt
fi
