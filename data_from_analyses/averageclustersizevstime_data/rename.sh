#!/bin/bash

#filenames=`find *dt0.002_combined*`
filenames=`find *gel_l860_ls1400_vfr*`

#for file in $filenames ; do mv $file ${file//dt0.002_combined/dt0.002_gs0.001_combined} ; done
for file in $filenames ; do mv $file ${file//gel_l860_ls1400_vfr/gel_l860_vfr} ; done
