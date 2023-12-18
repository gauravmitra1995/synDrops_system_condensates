This folder contains the code for dynamic binding / unbinding between particles. 

**(Singularity files have to be present inside the folder *dybond/singularity/* for the compilation of the plugin to work properly!!!!)**

**Download the singularity files from here:**

https://drive.google.com/file/d/15vH54mzLhiscVWF9_P6IRNhIRQntYn3Q/view?usp=sharing


**After downloading the singularity files (present as a zip file in the link above), unzip them and save the files individually in the folder *dybond/singularity/* before running anything else.**


**Information about important files in this folder**

The folder *nvidia-hoomd-2.9.6* contains the script *DyBondUpdater.cc* which basically contains the code for this dynamic binding. It also contains a header file *DyBondUpdater.h* where all the important global variables, functions and constructors are defined (which are used in the .cc script). A python script *update.py* is also present which activates the c++ DyBondUpdater from cpp module by initializing the dybond plugin and also feeds parameters to the DyBondUpdater written in C++ which are sent from the python simulation run script. 


**Instructions for using the dybond plugin (with singularity present)**

The dynamic bonding plugin has already been compiled and the *.so* file has been generated inside the folder *nvidia-hoomd-2.9.6*. 
The wrapper script to finally use for running simulations or performing analyses is: *run-hoomd2.9.6.bash* present inside the folder *dybond*.



