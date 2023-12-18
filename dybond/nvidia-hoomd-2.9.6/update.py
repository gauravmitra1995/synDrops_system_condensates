# Copyright (c) 2009-2016 The Regents of the University of Michigan
# This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

# this simple python interface just activates the c++ DyBondUpdater from cpp module

# Next, since we are extending an updater, we need to bring in the base class updater and some other parts from
# hoomd_script
import hoomd
from hoomd import _hoomd
from hoomd.md import _md

# First, we need to import the C++ module. It has the same name as this module (dybond_plugin) but with an underscore
# in front
from hoomd.dybond_plugin import _dybond_plugin

from hoomd.md import nlist as nl # to avoid naming conflicts
import copy;
# Pick a random particle and search for candidate neighbours to bond to and make the bonds if possible
#
# Every \a period time steps, the bonding attempt is made
#


class dybond(hoomd.update._updater):
    # Initialize the dybond plugin
    #
    # \param period: Dynamic bond creation will be attempted every \a period time steps
    #        nlist : The neighbourlist to be used by the updater.
    #        group : Group used to calculate the temperature.
    # \b Examples:
    # \code
    # dybond_plugin.update.dybond()
    # bondr = dybond_plugin.update.dybond(period=10)
    # \endcode
    #
    # \a period can be a function: see \ref variable_period_docs for details
    def __init__(self, nlist, group, period=1):
        '''
        Creates an instance of the dybond updater.
        :param nlist: The neighbourlist to be used by the updater.
        :param group: Group used to calculate the temperature.
        :param period: Number of time steps between two bonding attempts
        NOTE: Not tested for running using multiple GPU's using MPI
        '''
        hoomd.util.print_status_line();

        # initialize base class
        hoomd.update._updater.__init__(self);

        #store the neighborlist for later use
        self.bond_type_dict = []
        self.nlist = nlist
        #self.nlist.subscribe(lambda:self.get_rcut())

        #store the group information for later use
        self.group = group

        # create the compute thermo
        
        if group is hoomd.context.current.group_all:
            group_copy = copy.copy(group);
            group_copy.name = "__nvt_all";
            hoomd.util.quiet_status();
            self.thermo = hoomd.compute.thermo(group_copy);
            self.thermo.cpp_compute.setLoggingEnabled(False);
            hoomd.util.unquiet_status();
        else:
            self.thermo = hoomd.compute._get_unique_thermo(group=group);

        # initialize the reflected c++ class
        #if not hoomd.context.exec_conf.isCUDAEnabled():
        self.cpp_updater = _dybond_plugin.DyBondUpdater(hoomd.context.current.system_definition,
                                                        self.nlist.cpp_nlist,
                                                        self.thermo.cpp_compute);
        #else:
        #    self.cpp_updater = _dybond_plugin.DyBondUpdaterGPU(hoomd.context.current.system_definition, self.nlist.cpp_nlist);
        self.setupUpdater(period);

    def set_params(self,
                   bond_type,
                   A,B,nondybonds,r0,
                   rmin,rmax,kspring,flag_temperature,flag_force,
                   flag_metropolis,Tmelt,alpha,kon_init,kon_melt,koff_init,koff_melt,checksteps,dt,particle_clusterids,self_avoiding,userseed):
        '''
        Sets parameters to the Dybond updater. These parameters are in principle specific to each bond type.
        Parameters general for all bond types are set in the initializer.
        :param bond_type: Name of bond type e.g. 'C-D' or 'polymer'
        :param A: Name of the first particle type participating in the bond e.g. 'C'
        :param B: Name of the second particle type participating in the bond e.g. 'D'
        :return:
        '''
        self.cpp_updater.set_params(bond_type, A, B,
                                    nondybonds,r0,rmin,rmax,kspring,flag_temperature,flag_force,flag_metropolis,Tmelt,alpha,kon_init,kon_melt,koff_init,koff_melt,checksteps,
                                    dt,particle_clusterids,self_avoiding,userseed)
     
        assert (kon_init*checksteps*dt <= 1.0),"Error: kon_init is too high; P_on (initial) should be <= 1.0"
        assert (kon_melt == 0),"Error: kon_melt should be = 0"
        assert (kon_init+kon_melt-koff_init == koff_melt),"Error: koff_melt not calculated correctly; koff_melt = kon_init + kon_melt - koff_init"
        typeA = hoomd.context.current.system_definition.getParticleData().getTypeByName(A)
        typeB = hoomd.context.current.system_definition.getParticleData().getTypeByName(B)
        
        self.bond_type_dict.append([bond_type,typeA,typeB,nondybonds,r0,rmin,rmax,kspring,flag_temperature,flag_force,flag_metropolis,Tmelt,alpha,kon_init,kon_melt,koff_init,koff_melt,checksteps,dt,particle_clusterids,self_avoiding,userseed])
        # update the neighbor list for this bond type
        #self.nlist.update_rcut()

    ## \internal
    # \brief Get the r_cut pair dictionary
    # \returns The rcut(i,j) dict if atleast one bond type is defined and None otherwise
    def get_rcut(self):
        if len(self.bond_type_dict) < 1:
            hoomd.context.msg.warning("No cut offs are specified for neighborlist used in DyBondPlugin\n")
            return None
        # go through the list of bond types scheduled for dynamic bonding
        r_cut_dict = nl.rcut();
        for bond_data in self.bond_type_dict:
            # set the r_cut value
            r_cut_dict.set_pair(bond_data[1],bond_data[2],bond_data[3]);

        return r_cut_dict;
