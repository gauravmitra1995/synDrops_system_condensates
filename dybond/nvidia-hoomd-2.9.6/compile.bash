#!/bin/bash

/usr/bin/c++ -DCUDA_ARCH=35 -DCUSOLVER_AVAILABLE -DEIGEN_MPL2_ONLY -DENABLE_CUDA -DENABLE_HPMC_MIXED_PRECISION -DSINGLE_PRECISION -D_REENTRANT -D_dybond_plugin_EXPORTS -DNO_IMPORT_ARRAY -I/usr/local/cuda/include -I/ext3/miniconda3/include/python3.8 -I/ext3/hoomd-2.9.6-install/hoomd-v2.9.6 -I/ext3/hoomd-2.9.6-install/hoomd-v2.9.6/hoomd/extern/cereal/include -I/ext3/hoomd-2.9.6-install/build/hoomd/include  -Wall -Wno-unknown-pragmas -Wno-deprecated-declarations -O3 -DNDEBUG -fPIC   -std=gnu++11 -o DyBondUpdater.cc.o -c DyBondUpdater.cc

/usr/bin/c++ -fPIC  -Wall -Wno-unknown-pragmas -Wno-deprecated-declarations -O3 -DNDEBUG  -shared -Wl,-soname,_dybond_plugin.so \
	     -o _dybond_plugin.so \
	     DyBondUpdater.cc.o \
	     /ext3/hoomd-2.9.6-install/build/hoomd/dybond_plugin/CMakeFiles/_dybond_plugin.dir/module.cc.o \
	     /ext3/hoomd-2.9.6-install/build/hoomd/dybond_plugin/CMakeFiles/cuda_compile_1.dir/cuda_compile_1_generated_DyBondUpdater.cu.o \
	     -Wl,-rpath,/usr/local/cuda/lib64:/ext3/hoomd-2.9.6-install/build/hoomd/md:/ext3/hoomd-2.9.6-install/build/hoomd/extern/neighbor/neighbor:/ext3/hoomd-2.9.6-install/build/hoomd: /usr/local/cuda/lib64/libcudart_static.a -lpthread -ldl -lrt /usr/local/cuda/lib64/libcufft.so /usr/local/cuda/lib64/libcurand.so -ldl -lpthread \
	     /ext3/hoomd/hoomd/md/_md.cpython-38-x86_64-linux-gnu.so \
	     -lrt /usr/local/cuda/lib64/libcufft.so /usr/local/cuda/lib64/libcurand.so \
	     /ext3/hoomd/hoomd/_hoomd.cpython-38-x86_64-linux-gnu.so \
	     /ext3/hoomd/hoomd/libquickhull.so \
	     -Wl,-rpath-link,/ext3/hoomd-2.9.6-install/build/hoomd/extern/neighbor/neighbor 
