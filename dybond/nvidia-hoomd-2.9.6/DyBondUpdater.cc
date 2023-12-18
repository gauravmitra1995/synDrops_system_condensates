#define USE_DIST_SORTED_NLIST
#define USE_MERSENNE_TWISTER_ENGINE
#include "DyBondUpdater.h"
#ifdef ENABLE_CUDA
#include "DyBondUpdater.cuh"
#endif
#include "hoomd/GPUArray.h"
#include <math.h>
#include <ctime>
#include <bits/stdc++.h> 
#include <iostream>
#include <string>
#include <random>

#define random() 1.0*rand()/RAND_MAX

// windows defines a macro min and max
#undef min
#undef max

using namespace hoomd;
using namespace std;

// start the random number generator
std::mt19937_64 rng(unsigned(time(NULL)));
std::uniform_real_distribution<double> unif(0, 1);

/*
    Dynamic Bond Updater creates bonds between two types of particles.
*/


/* \param sysdef: System to dynamically create bonds
    \param nlist : neighbour list to use for finding particles to bond to
    \param thermo: ComputeThermo object to get temperature of the system
    \param group : The group of particles this updater method is to work on
*/
DyBondUpdater::DyBondUpdater(std::shared_ptr<SystemDefinition> sysdef, 
                             std::shared_ptr<NeighborList> nlist,std::shared_ptr<ComputeThermo> thermo)
        : Updater(sysdef), m_nlist(nlist), m_thermo(thermo)
    {
    // access the bond data for later use
    m_bond_data = m_sysdef->getBondData();
    assert(m_nlist);
    }

/* Defines a bond type between two particle types that should be created by the updater
*/
void DyBondUpdater::set_params(std::string bond_type,
                               std::string A,std::string B,Scalar nondybonds,Scalar r0,
                               Scalar rmin,Scalar rmax,Scalar kspring,unsigned int flag_temperature, double flag_force, 
                               unsigned int flag_metropolis,Scalar Tmelt,Scalar alpha,Scalar kon_init,Scalar kon_melt,Scalar koff_init,Scalar koff_melt,
                        Scalar checksteps,Scalar dt,std::string particle_clusterids,unsigned int self_avoiding, unsigned int userseed)
    {
    if (m_bond_type.size() < 1)
        {	
        m_seed=userseed;
        rng.seed(m_seed);

	first_type=A;
	second_type=B;
	num_nondybonds=nondybonds;

	b_type = m_bond_data->getTypeByName(bond_type);
	A_type = m_pdata->getTypeByName(A);
	B_type = m_pdata->getTypeByName(B);

	cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
	cout<<"Using DynBondUpdater latest updated version (October 2023) after resolving the most recent bug with Pon (with metropolis)"<<endl;
        cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;

	cout<<"First bondable particle type: "<<first_type.c_str()<<endl;
	cout<<"Second bondable particle type: "<<second_type.c_str()<<endl;
	
	m_rmin=rmin;
	m_rmax=rmax;
	m_r0=r0;
	m_kspring=kspring;
	m_checksteps=checksteps;
	m_Tmelt=Tmelt;
	m_alpha=alpha;
	m_dt=dt;
        m_koninit=kon_init;
        m_konmelt=kon_melt;
        m_koffinit=koff_init;
        m_koffmelt=koff_melt;
        m_particle_clusterids=particle_clusterids;
        flag_T=flag_temperature;
        flag_F=flag_force;
        flag_MP=flag_metropolis;
        flag_selfavoiding=self_avoiding;

       
        std::string fragment = ""; 
        for (auto x : m_particle_clusterids) 
        { 
           if (x == ' ') 
           { 
              clusterids.push_back(std::stoi(fragment));
              fragment = ""; 
           } 
           else
           { 
              fragment = fragment + x; 
           } 
        }

        clusterids.push_back(std::stoi(fragment)); 

        if(A_type==B_type)
        {
           cout<<"Particle types for bonding are similar (homogeneous)"<<endl;
        }
        else
        {
           cout<<"Particle types for bonding are not same at all"<<endl;
        }
        
        if (m_bond_type.find(b_type) == m_bond_type.end())
            m_bond_type[b_type] = std::make_pair(A_type,B_type);
        else
            m_exec_conf->msg->warning() << "Bond type: " << b_type << "is already defined!";
            
        
        if (m_bonds_made.find(b_type) == m_bonds_made.end())
            m_bonds_made[b_type] = 0; 
        
	//resetting bonds made count in the beginning.
        //populate the pair that stores the bond type and their participants as <bond_type,participant1 type, participant2 type>


        init_dictionaries_from_system(b_type);
        }
        else
        {
        m_exec_conf->msg->warning() << "The set_params function of dybond updater was used more than once."
                                    << " As of now, the dybond updater is only tested for one bond type."
                                    << " Hence only the first call to set_params has been used. Others were ignored\n";
        }
    }

/* \param bond_type: bond type to search from system and populate dictionaries
    This is done so that if the initial condition contains some bonds of type that
    are being made by the DyBondUpdater, it needs to know those bonds so that we
    do not make duplicate bonds and the logged quantities reflect the actual numbers.
    This is required because we use two dictionaries to speed up look up during the actual simulation.
    They are
    1) m_bond_data = { bond type : number of bonds of this type }
    2) m_rank_dict = { particle id : number of bonds made by this particle }
*/


void DyBondUpdater::init_dictionaries_from_system(unsigned int bond_type)
    {

    ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_tag(m_pdata->getTags(), access_location::host, access_mode::read);
    ArrayHandle< unsigned int > h_rtag(m_pdata->getRTags(), access_location::host, access_mode::read);
 
    ArrayHandle<typename BondData::members_t> h_bonds(m_bond_data->getMembersArray(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int>  h_bond_tags(m_bond_data->getTags(), access_location::host, access_mode::read);

    initialbonds=num_nondybonds;
    cout<<"No of initial non-dynamic bonds:"<<num_nondybonds<<endl;
    cout<<"Given Dynamic Bond type:"<<bond_type<<endl;


    unsigned int num_bonds = 0;
    // for each of the bonds
    const unsigned int size = (unsigned int)m_bond_data->getN();
    cout<<"No of actual initial bonds:"<<size<<endl;

    total_particles = (unsigned int) m_pdata->getN();
    cout<<"Total number of particles:"<<total_particles<<endl;

    for (auto& current_bond_type_entry : m_bond_type)
    {
       auto curr_b_type = current_bond_type_entry.first;
       cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
       cout<<"curr_bond_type:"<<curr_b_type<<endl;

       for (unsigned int i=0;i<total_particles;i++)
       {

       m_rank_dict[i] = 0;
       m_tbind[i] = -1;
       m_tunbind[i] = -1;
   
       unsigned int p_from_type;
       unsigned int p_from_idx=i;
       unsigned int p_from_tag=h_rtag.data[p_from_idx];

       p_from_type = __scalar_as_int(h_pos.data[p_from_tag].w);

       if(!is_bondable_particle(curr_b_type, p_from_type))
       {
         continue;
       }
       else
       {
       dynbondparticles.push_back(i);   //Creating a vector consisting of only the particles to be involved in dynamic bonding for the current bond type
       }

       }

       cout<<"No of eligible particles for bonding:"<<dynbondparticles.size()<<endl;

    }

    for (unsigned int i = 0; i < size; i++)
        {
        auto curr_b_type = m_bond_data->getTypeByIndex(i);
	
        if(curr_b_type!=bond_type)
		continue;
        // lookup the tag of each of the particles participating in the bond
        const BondData::members_t bond = m_bond_data->getMembersByIndex(i);
        assert(bond.tag[0] < m_pdata->getN());
        assert(bond.tag[1] < m_pdata->getN());

        auto a_bond_rank=1;
        auto b_bond_rank=1;
        
        m_rank_dict[bond.tag[0]]=a_bond_rank;
        m_rank_dict[bond.tag[1]]=b_bond_rank;
	
        num_bonds++;
        }
	
    m_bonds_made[bond_type] = num_bonds;       

    cout<<"No of bonds of type "<<bond_type<<" is "<<m_bonds_made[bond_type]<<endl;
   }
    


/* Returns the number of times this particle has been bonded already
*/

unsigned int DyBondUpdater::bond_rank(unsigned int p_idx)
    {
    if (m_rank_dict.find(p_idx) == m_rank_dict.end())
    {
        m_rank_dict[p_idx]=0;
    }
    return m_rank_dict[p_idx];
    }


/* Checks if particle type "p_type" is a participating type in bond type "b_type"
 */

bool DyBondUpdater::is_bondable_particle(unsigned int b_type, unsigned int p_type)
    { 
    return (m_bond_type[b_type].first == p_type || m_bond_type[b_type].second == p_type);       
    }

//Calculate koff (no force dependence)
long double DyBondUpdater::calc_koff_noforce(double m_Tmelt,double m_alpha,double curr_T,double koff_init,double koff_melt,unsigned int flag_T)
     {
        long double g=0.5*(tanh(m_alpha*(curr_T-m_Tmelt))+1.0);
        long double koff=0.0;       
                         
        if(flag_T == 0)
        {
            koff=koff_init;
        }
        else if(flag_T == 1)
        {
            koff=koff_init*(1.0-g)+koff_melt*g;
        }

        return koff;
     }

//Calculate koff (with force dependence)
long double DyBondUpdater::calc_koff_force(double curr_T,double koff_init,double distance,double m_r0,double m_kspring,double flag_F)
     {
	
        long double koff=0.0;

        if(flag_F == 0.0)
        {
            koff=koff_init;
        }
        else      // Implement slip bond dependence 
        {
            double a=flag_F;
	    long double Force=m_kspring*(fabs(distance-m_r0));
            long double d=sqrt(curr_T/m_kspring);
            koff=koff_init*exp((a*Force*d)/(curr_T));
        }
        return koff;
     }

//Calculate kon
long double DyBondUpdater::calc_kon(double m_Tmelt,double m_alpha,double curr_T,double kon_init,double kon_melt,unsigned int flag_T)
     {
        long double f=0.5*(tanh(m_alpha*(curr_T-m_Tmelt))+1.0);        
        long double kon=0.0;

        if(flag_T == 0)
        {
            kon=kon_init;
        }
        else if(flag_T == 1)
        {
            kon=kon_init*(1.0-f)+kon_melt*f;
        }
	
        return kon;
     }

Scalar DyBondUpdater::calculate_distance(unsigned int p_from_tag,unsigned int p_to_tag)
     {
        // access the particle data
        ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(), access_location::host, access_mode::read);
        ArrayHandle< unsigned int > h_rtag(m_pdata->getRTags(), access_location::host, access_mode::read);
        ArrayHandle< unsigned int > h_tag(m_pdata->getTags(), access_location::host, access_mode::read);
        
        BoxDim box = m_pdata->getBox();

        vec3<Scalar> pi(h_pos.data[p_from_tag].x, h_pos.data[p_from_tag].y, h_pos.data[p_from_tag].z);
        vec3<Scalar> pj(h_pos.data[p_to_tag].x, h_pos.data[p_to_tag].y, h_pos.data[p_to_tag].z);
        vec3<Scalar> dxScalar(pj - pi);
        
        // apply periodic boundary conditions
        dxScalar = vec3<Scalar>(box.minImage(vec_to_scalar3(dxScalar)));
        Scalar dist = sqrt(dxScalar.x*dxScalar.x + dxScalar.y*dxScalar.y + dxScalar.z*dxScalar.z);
        return dist;
     }

/* The main function in the updater code that gets called every "checksteps" period. 
 *  Computes the neighbourlist, temperature and uses particle information to make bonds with a probability.
    \param timestep Current time step of the simulation
*/

void DyBondUpdater::update(unsigned int timestep)
    {
        if (m_prof) m_prof->push("DyBondUpdater");
        // start by updating the neighborlisttest_dybond.py
        m_nlist->compute(timestep);
        // compute the current thermodynamic properties and get the temperature
        m_thermo->compute(timestep);
        curr_T = m_thermo->getTranslationalTemperature();

        // access the neighbor list
        ArrayHandle<unsigned int> h_n_neigh(m_nlist->getNNeighArray(), access_location::host, access_mode::read);
        ArrayHandle<unsigned int> h_nlist(m_nlist->getNListArray(), access_location::host, access_mode::read);
        ArrayHandle<unsigned int> h_head_list(m_nlist->getHeadList(), access_location::host, access_mode::read);
        // access the particle data
        ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(), access_location::host, access_mode::read);
        ArrayHandle< unsigned int > h_rtag(m_pdata->getRTags(), access_location::host, access_mode::read);
        ArrayHandle< unsigned int > h_tag(m_pdata->getTags(), access_location::host, access_mode::read);

        //sanity check
        assert(m_bond_type.size() > 0);
	
        //try bonding and unbonding as requested by the user
        for (auto& current_bond_type_entry : m_bond_type)
        {

               // Access the CPU bond table for reading
               ArrayHandle<typename BondData::members_t> h_bonds(m_bond_data->getMembersArray(), access_location::host, access_mode::read);
               ArrayHandle<unsigned int>  h_bond_tags(m_bond_data->getTags(), access_location::host, access_mode::read);

               auto curr_b_type = current_bond_type_entry.first;
               auto curr_b_pair = current_bond_type_entry.second; 

               //UNBINDING:

               const unsigned int numberofbonds = (unsigned int)m_bond_data->getN();
	       bonds_to_remove.clear();
 
               for (unsigned int bond_number = initialbonds; bond_number < numberofbonds; bond_number++) 
               {

		 // look up the tag of each of the particles participating in the bond
		 const BondData::members_t bond = m_bond_data->getMembersByIndex(bond_number);
		 assert(bond.tag[0] < m_pdata->getN());
		 assert(bond.tag[1] < m_pdata->getN());          

		 unsigned int tag_first=h_rtag.data[bond.tag[0]];
		 unsigned int tag_second=h_rtag.data[bond.tag[1]];

		 assert(tag_first <= m_pdata->getMaximumTag());
		 assert(tag_second <= m_pdata->getMaximumTag());

                 Scalar distance=calculate_distance(tag_first,tag_second);

		 unsigned int type_first=__scalar_as_int(h_pos.data[tag_first].w);
		 unsigned int type_second=__scalar_as_int(h_pos.data[tag_second].w);

		 unsigned int bondrank_p1=bond_rank(bond.tag[0]);
		 unsigned int bondrank_p2=bond_rank(bond.tag[1]);
                 
		 //Checking if the two particles we are trying to unbind actually belong to the dynamic bond type under consideration or not. 
		 if(!((type_first==curr_b_pair.first && type_second==curr_b_pair.second) || (type_second==curr_b_pair.first && type_first==curr_b_pair.second)))
	         {
			 continue;
		 }

                 m_Poff=0.0;

		 if(flag_F==0)
		 {
		 m_koff=calc_koff_noforce(m_Tmelt,m_alpha,curr_T,m_koffinit,m_koffmelt,flag_T);
		 }
		 else
		 {
		 m_koff=calc_koff_force(curr_T,m_koffinit,distance,m_r0,m_kspring,flag_F);
		 }
		 m_Poff=m_koff*m_checksteps*m_dt;      //Calculating Poff from koff

		 double r;
		 r=unif(rng);
                 
		 if(r<m_Poff)
	         {
		      cout<<"Unbonding particles "<<bond.tag[0]<<" (bond rank "<<bondrank_p1<<")"<<" and "<<bond.tag[1]<<" (bond rank "<<bondrank_p2<<")"<<" of ptypes "<<type_first<<" and "<<type_second<<" distant "<<distance<<" with T:"<<curr_T<<", Poff:"<<m_Poff<<" and random number r:"<<r<<endl;
		      m_rank_dict[bond.tag[0]]=0;
		      m_rank_dict[bond.tag[1]]=0;
		      m_bonds_made[curr_b_type]--;       
		      m_tunbind[bond.tag[0]]=timestep;
		      m_tunbind[bond.tag[1]]=timestep;
		      bonds_to_remove.push_back(bond_number);
		      
	         }
		 
	       }
              
              //remove bonds in reverse order
              for ( int i=bonds_to_remove.size(); i>0; i-- ){
                  m_bond_data->removeBondedGroup(h_bond_tags.data[bonds_to_remove[i-1]]);
              }

              
	      
	      //BINDING:
	       
              bonds_to_add_dist_dict.clear();
              bonds_to_add_idx_dict.clear();
              bonds_to_add_to_idx_dict.clear();
              bonds_to_add_idx_list.clear();

               for (unsigned int i=0;i<dynbondparticles.size();i++) //loop over all eligible particles to bind 
               {
		   
                  unsigned int p_to_type;
                  unsigned int p_from_type;
                  unsigned int p_from_idx=dynbondparticles[i];    

                  unsigned int p_from_tag=h_rtag.data[p_from_idx];
                  p_from_type = __scalar_as_int(h_pos.data[p_from_tag].w); 
                  unsigned int state1 = bond_rank(p_from_idx);
                  
		  //if this particle 'p_from_idx' is already in a bond proposal, skip:

                  if (bonds_to_add_to_idx_dict.find(p_from_idx) != bonds_to_add_to_idx_dict.end()){
			  continue;
                  }

		  if(bonds_to_add_idx_dict.find(p_from_idx) != bonds_to_add_idx_dict.end()){
			  continue;

	          }
		  	  		  
                  m_Pon=0.0;

                  if(state1==0) //if bond rank of particle is 0, then attempt binding 
                  {

                    if(m_tunbind[p_from_idx]==timestep) //Skip this particle since it was already unbound at the same time
                    {
                    continue;
                    }

                    //pick the complementary particle to bond to
                    auto p1_type = std::get<0>(curr_b_pair);
                    auto p2_type = std::get<1>(curr_b_pair);
                    p_to_type = (p_from_type==p1_type) ? p2_type : p1_type;

                    // sanity check
                    assert(p_from_type < m_pdata->getNTypes());

                    // loop over all of the neighbors of this particle and try bonding
                    const unsigned int size = (unsigned int)h_n_neigh.data[p_from_tag];
                    const unsigned int head_i = h_head_list.data[p_from_tag];
                    std::vector<Scalar> distances;
		    std::vector<Scalar> distances_minus_r0;
                    std::vector<unsigned int> p_to_idxs;

                    if(size!=0)
                    {
	              
                      for (unsigned int j = 0; j < size; j++) //over all neighbors for the particle
                      {
		       
                       unsigned int p_to_tag = h_nlist.data[head_i + j];
                       auto p_to_idx = h_tag.data[p_to_tag];
		       Scalar dist = calculate_distance(p_from_tag,p_to_tag);
                       
                       auto state2 = bond_rank(p_to_idx);
                       unsigned int p_temp_type;
                       p_temp_type = __scalar_as_int(h_pos.data[p_to_tag].w);
                       auto bondable_type = (p_temp_type == p_to_type);
                       auto within_cutoff = (dist>m_rmin && dist<m_rmax);
                       
		       signed int cluster_p1=clusterids[p_from_idx];
		       signed int cluster_p2=clusterids[p_to_idx];
		       signed int difference=cluster_p1-cluster_p2;
                      

		       if(state2!=0)  //check if the particle to be bonded to is itself already bonded or not
		       {
			   continue;
		       }
		       if(!(bondable_type)) //if particle to be bonded to belongs to the complementary type
		       {
			   continue;
		       }
		       if(clusterids[p_from_idx]==clusterids[p_to_idx]) //if they belong to the same cluster
		       {
			   continue;
		       }
		       else
		       {
                          if(flag_selfavoiding==1 and abs(difference)!=1) //check if the difference in cluster indices obey 'self-avoiding' criterion 
		          {
		             continue;
	                  }
		       }

                       if (within_cutoff)  //if a particle falls within the bonding range, then add those particle indices and distances to vector(s) 
                       {
                          distances.push_back(dist);
			  distances_minus_r0.push_back(fabs(dist-m_r0));
                          p_to_idxs.push_back(p_to_idx);
                       }
                      

                      }
		       
                      if(distances_minus_r0.size()!=0)
                      {
			 //Choose the best neighbor for binding based on minimum of fabs(dist-r0):
			 	
			 unsigned int min_index = min_element(distances_minus_r0.begin(), distances_minus_r0.end())-distances_minus_r0.begin();
			 unsigned int closest_particle_idx = p_to_idxs[min_index];                        

                                                    
			  //if this particle that p_from_idx is to be bonded to is already in a bond proposal, skip:

			  if (bonds_to_add_to_idx_dict.find(closest_particle_idx) != bonds_to_add_to_idx_dict.end()){
				  continue;
			  }

			  if(bonds_to_add_idx_dict.find(closest_particle_idx) != bonds_to_add_idx_dict.end()){
				  continue;

			  }
                         
                         //add this pair to the binding list
			 
                         bonds_to_add_idx_dict[p_from_idx]=closest_particle_idx;
                         bonds_to_add_to_idx_dict[closest_particle_idx]=p_from_idx;
                         bonds_to_add_dist_dict[p_from_idx]=distances[min_index];
                         bonds_to_add_idx_list.push_back(p_from_idx);
		      }
                    }

                  }
	       }		 

               m_kon=calc_kon(m_Tmelt,m_alpha,curr_T,m_koninit,m_konmelt,flag_T);
               m_Pon_default=m_kon*m_checksteps*m_dt;                  //Calculating Pon from kon
                 
	     
	       for(unsigned int k=0;k<bonds_to_add_idx_list.size();k++)
	       {

		 unsigned int p_from_idx = bonds_to_add_idx_list[k];
		 unsigned int p_to_idx = bonds_to_add_idx_dict[p_from_idx];
		 double distance = bonds_to_add_dist_dict[p_from_idx];

		 unsigned int p_from_tag=h_rtag.data[p_from_idx];
                 unsigned int p_from_type = __scalar_as_int(h_pos.data[p_from_tag].w);

		 unsigned int p_to_tag=h_rtag.data[p_to_idx];
		 unsigned int p_to_type = __scalar_as_int(h_pos.data[p_to_tag].w);

		 unsigned int bondrank_p1=bond_rank(p_from_idx);
	         unsigned int bondrank_p2=bond_rank(p_to_idx);

		 //add Metropolis criterion to satisfy detailed balance

		 m_Pon=(m_Pon_default)*exp(-0.5*flag_MP*m_kspring/curr_T*pow(distance-m_r0,2)); //resolving the earlier bug here (by introducing a default Pon variable)  Pon would be changing everytime inside this loop
		 
		 double r;
		 r=unif(rng);
                 
		 //Print warning 
		 if(bondrank_p1==1 or bondrank_p2==1)
		 {
			 continue;
	         }
	         
		 if(r<m_Pon)
		 {
		  m_tbind[p_from_idx]=timestep;
		  m_tbind[p_to_idx]=timestep;

		  cout<<"Bonding particles "<<p_from_idx<<" (bond rank "<<bondrank_p1<<")"<<" and "<<p_to_idx<<" (bond rank "<<bondrank_p2<<")"<<" of ptypes "<<p_from_type<<" and "<<p_to_type<<" distant "<<distance<<" with T:"<<curr_T<<", Pon:"<<m_Pon<<" and random number r:"<<r<<endl;
		  
		  m_bond_data->addBondedGroup(Bond(curr_b_type,p_from_idx,p_to_idx));
		  m_rank_dict[p_from_idx]=1;
		  m_rank_dict[p_to_idx]=1;
		  m_bonds_made[curr_b_type]++;
		 }       
         
               }
     
      }
     if(m_prof) m_prof->pop();
     }

std::vector< std::string > DyBondUpdater::getProvidedLogQuantities()
    {
    std::vector<std::string> ret;
    for(auto& btype : m_bond_type)
        {
        auto bond_type_str = m_bond_data->getNameByType(btype.first);
        }
    m_loggable_quantities = ret; //This is a deep copy. Any more loggable quantitites need to be appended.
    return ret;
    }

Scalar DyBondUpdater::getLogValue(const std::string& quantity, unsigned int timestep)
    {
    if (std::find(m_loggable_quantities.begin(),m_loggable_quantities.end(),quantity) != m_loggable_quantities.end())
        {
        unsigned int first = quantity.find('(')+1;
        unsigned int last = quantity.find_last_of(')');
        std::string this_b_name = quantity.substr(first,last-first);
        std::string this_quantity_name = quantity.substr(0,first-1);
        }
    else
        {
        m_exec_conf->msg->error() << "update.dybond updater: " << quantity
                                  << " is not a valid log quantity." << std::endl;
        throw std::runtime_error("Error getting log value");
        }
    return Scalar(0.0);
    }

/* Export the CPU updater to be visible in the python module
 */
void export_DyBondUpdater(pybind11::module& m)
    {
    pybind11::class_<DyBondUpdater, std::shared_ptr<DyBondUpdater> >(m, "DyBondUpdater", pybind11::base<Updater>())
        .def(pybind11::init<std::shared_ptr<SystemDefinition>,
        std::shared_ptr<NeighborList>,
        std::shared_ptr<ComputeThermo> >())
        .def("set_params", &DyBondUpdater::set_params)
    ;
    }	
