// Copyright (c) 2009-2016 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// inclusion guard
#ifndef _DYBOND_UPDATER_H_
#define _DYBOND_UPDATER_H_

/*! \file DyBondUpdater.h
*/
#include "hoomd/md/NeighborList.h"
#include "hoomd/ComputeThermo.h"
#include "hoomd/ParticleGroup.h"
#include <hoomd/Updater.h>
#include <map>
#include <random>
//#ifndef NVCC
#if !defined NVCC && !defined __HIP_PLATFORM_HCC__
#include <hoomd/extern/pybind/include/pybind11/pybind11.h>
#endif

//! An updater for dynamically creating and removing bonds between pairs of eligible particle types
class DyBondUpdater : public Updater
    {
    public:
        //! Constructor
        DyBondUpdater(std::shared_ptr<SystemDefinition> sysdef, std::shared_ptr<NeighborList> nlist,std::shared_ptr<ComputeThermo> thermo);
        //! Take one timestep forward
        virtual void update(unsigned int timestep);
        //! Returns a list of log quantities this compute calculates
        virtual std::vector< std::string > getProvidedLogQuantities(void);
        //! Calculates the requested log value and returns it
        virtual Scalar getLogValue(const std::string& quantity, unsigned int timestep);
        //! set the parameters that are not set in the constructor
        void set_params(std::string bond_type,
                        std::string A,
                        std::string B,
                        Scalar nondybonds,Scalar r0,
                        Scalar rmin,Scalar rmax,Scalar kspring,unsigned int flag_temperature,double flag_force,
                        unsigned int flag_metropolis,Scalar Tmelt,Scalar alpha,Scalar kon_init,Scalar kon_melt,Scalar koff_init,Scalar koff_melt,
                        Scalar checksteps,Scalar dt,std:: string particle_clusterids,unsigned int self_avoiding,unsigned int userseed);
        //std::map<unsigned int;

    protected:
        unsigned int m_seed;                                //user input seed to use
        unsigned int bond_rank(unsigned int p_idx);
        void init_dictionaries_from_system(unsigned int bond_type);
        bool is_bondable_particle(unsigned int b_type, unsigned int p_type);
        Scalar calculate_distance(unsigned int p_from_tag,unsigned int p_to_tag);
        long double calc_koff_noforce(double m_Tmelt,double m_alpha,double curr_T,double koff_init,double koff_melt,unsigned int flag_T);
	long double calc_koff_force(double curr_T,double koff_init,double distance,double m_r0,double m_kspring,double flag_F);
        long double calc_kon(double m_Tmelt,double m_alpha,double curr_T,double kon_init,double kon_melt,unsigned int flag_T);
        std::shared_ptr<BondData> m_bond_data;    //!< Bond data to use in computing bonds
        std::map<unsigned int,std::pair<unsigned int, unsigned int> > m_bond_type; //!< Dictionary of bond type and the participating particle types
        std::map<unsigned int, unsigned int> m_bonds_made; //!< Key is bond type and value is number of bonds made by the system.
        std::vector<std::string> m_loggable_quantities;
        std::shared_ptr<NeighborList> m_nlist;    //!< The neighborlist to use for the computation
        const std::shared_ptr<ComputeThermo> m_thermo; //!< compute for thermodynamic quantities
        Scalar m_rmin;
        Scalar m_rmax;               
        Scalar m_Pon_default; //default probability for binding
        Scalar m_Pon;//probability for binding after metropolis is applied
        Scalar m_Poff;//probability for unbinding
        Scalar m_kon;
        Scalar m_koff;
        Scalar m_r0;
        Scalar m_kspring;
        Scalar m_Tmelt;
        unsigned int flag_T;
        double flag_F;
	unsigned int flag_MP;
        unsigned int flag_selfavoiding;
        Scalar m_alpha;
        Scalar m_checksteps;
        Scalar m_dt;
        Scalar num_nondybonds;
        Scalar curr_T;
        Scalar m_koninit;
        Scalar m_konmelt;
        Scalar m_koffinit;
        Scalar m_koffmelt;
        std::vector<unsigned int> clusterids;
        std::vector<unsigned int> dynbondparticles;
	std::vector<unsigned int> bonds_to_remove;

        std::vector<unsigned int> bonds_to_add_idx_list;
        std::map<unsigned int, double> bonds_to_add_dist_dict;
        std::map<unsigned int, unsigned int> bonds_to_add_idx_dict;
        std::map<unsigned int, unsigned int> bonds_to_add_to_idx_dict;

        std::string first_type;
        std::string second_type;
        std::string m_particle_clusterids;
        std::map<unsigned int, unsigned int> m_rank_dict; //!< stores bond rank for particle ids to speed up lookup:
        std::map<unsigned int, float>m_tbind;
        std::map<unsigned int, float>m_tunbind;
        unsigned int A_type;
        unsigned int B_type;
        unsigned int b_type;
        unsigned int initialbonds;
        unsigned int total_particles;
    };

//! Export the DyBondUpdater class to python
void export_DyBondUpdater(pybind11::module& m);

#endif // _DYBOND_UPDATER_H_
