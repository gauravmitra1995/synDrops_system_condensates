import hoomd
import numpy as np
import os
import sys
import argparse
import hoomd.dybond_plugin as db
from hoomd import md
from potentials import softrepulsion

def gen_lattice(number_density,num_atoms):
    box_length = (number_density/num_atoms)**(-1./3)
    repeats = int(np.ceil((num_atoms/2.)**(1./3)))
    num_positions = 2*(repeats**3)

    lattice_spacing = (number_density/num_positions*num_atoms)**(-1./3)
    #lattice_spacing = (number_density/2)**(-1./3)

    from ase.lattice.cubic import SimpleCubicFactory
    class CsClFactory(SimpleCubicFactory):
        "A factory for creating CsCl"
        bravais_basis = [[0, 0, 0], [0.5, 0.5, 0.5], ]
        element_basis = (0,1)
    
    fact = CsCl = Rocksalt = CsClFactory()
    
    
    atoms = fact(directions=[[1,0,0], [0,1,0], [0,0,1]],
                size=(repeats,repeats,repeats), symbol=['N','P'],
                latticeconstant=lattice_spacing)
    positions = atoms.get_positions()
    np.random.shuffle(positions)
    positions -= box_length/2. * np.ones(positions.shape)
    return positions


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--gpu",default=False,action="store_true")
    parser.add_argument("--outprefix",help="Start of output simulation name (default: %(default)s)",default=None,type=str,required=True)
    parser.add_argument("--nsteps",help="Number of steps of simulation to run (default: %(default)s)",default=None,type=int,required=True)
    parser.add_argument("--dt",help="Simulation time step (default: %(default)s)",default=0.001,type=float)
    parser.add_argument("--koff",help="Off rate for unbinding in time units (default: %(default)s)",default=0.001,type=float)
    parser.add_argument("--sphere_repulsion",help="Repulsion of spheres including crowders in units of kT, for soft potential (default: %(default)s)",default=100,type=float)
    parser.add_argument("--rod_repulsion",help="Repulsion of spheres including crowders in units of kT, for soft potential (default: %(default)s)",default=500,type=float)
#box_length=860.
    parser.add_argument("--crowder_temperature",help="Temperature of crowders, relative to 1.0 (default: %(default)s)",default=1.0,type=float)
    parser.add_argument("--box_length",help="Side length of simulation box in nm (default: %(default)s)",default=860,type=float)
    parser.add_argument("--binding_distance",help="Maximum distance where binding can take place, in nm (default: %(default)s)",default=0.9,type=float)
    parser.add_argument("--number_rods",help="Number of rod proteins (default: %(default)s)",default=1170,type=int)
    parser.add_argument("--number_linkers",help="Number of rod proteins (default: %(default)s)",default=390,type=int)
    parser.add_argument("--volume_fraction_crowders",help="Fraction of volume taken up by crowders (default: %(default)s)",default=0,type=float)
    parser.add_argument("--seed",help="Random seed (default: %(default)s)",default=1,type=int)
    parser.add_argument("--dump_frequency",help="Dump frequency in time units(default: %(default)s)",type=int,default=20)
    parser.add_argument("--lj_epsilon",help="Epsilon for LJ attraction (default: %(default))",type=float,default=0)
    parser.add_argument("-i","--inputfile",default=None,help="Input file for continuing (gsd file, optional)",type=str)

    args = parser.parse_args()
    print("Command line arguments given:")
    print(args)

    locals().update(vars(args))
    if gpu:
        hoomd.context.initialize("--mode=gpu")
        print("Running on the GPU")
    else:
        hoomd.context.initialize("--mode=cpu")

    gsd_file = outprefix+".gsd"
    log_file = outprefix+".log"

    np.random.seed(seed)
    dump_frames = int(dump_frequency/dt)

    rod_segments = 2
    diameter_rod = 23.4/float(rod_segments)
    diameter_linker = 12.6
    diameter_crowder = 30.
    radius_sticker = 1
    bond_length=2*radius_sticker
    spring_constant=20

    if inputfile is None:
        volume_in_nm = (box_length)**3
        volume_in_microns = (box_length/1000.)**3
        
        #in nm 3
        volume_crowder = 4./3*np.pi*((diameter_crowder/2.)**3)
        number_crowders = int(volume_fraction_crowders*volume_in_nm/volume_crowder)
        number_particles = number_rods + number_linkers + number_crowders
        print("Starting the system with %i linkers, %i rods, %i crowders"%(number_linkers, number_rods, number_crowders))
    
        particle_types = particle_types=['rod', 'linker','crowder','r','l','s']
        snapshot = hoomd.data.make_snapshot(N=number_particles, box=hoomd.data.boxdim(L=box_length,dimensions=3), particle_types=particle_types, 
                      bond_types=['r-l'])
    
        number_density = number_particles / volume_in_nm
        #random_lattice_positions = gen_lattice(number_density,number_particles)
        random_lattice_positions = (np.random.random(size=(number_particles,3))-0.5)*box_length
        
        for i in range(number_particles):
            if i < number_rods:
                snapshot.particles.typeid[i]=0
                snapshot.particles.diameter[i]=diameter_rod
                snapshot.particles.moment_inertia[i] = np.array([1,1,0])*0.1
            elif i >= number_rods and i < number_rods + number_linkers:
                snapshot.particles.typeid[i]=1
                snapshot.particles.diameter[i]=diameter_linker
                snapshot.particles.moment_inertia[i] = np.array([1,1,1])*1
            else:
                snapshot.particles.typeid[i]=2
                snapshot.particles.diameter[i]=diameter_crowder
                snapshot.particles.moment_inertia[i] = np.array([1,1,1])
            snapshot.particles.mass[i]=1
            snapshot.particles.position[i] = random_lattice_positions[i]
        
        snapshot.bonds.resize(0)
        #snapshot.bonds.group[0]=[0,1]
        system=hoomd.init.read_snapshot(snapshot)
        
        # Add constituent particles of type rod and create the rods
        system.particles.types.add('R');
        
        rigid = hoomd.md.constrain.rigid();
        rod_positions = [0]
        for i in range(rod_segments+1):
            if i == 0 or i == rod_segments: 
                shift = radius_sticker
            else:
                shift = 0
            rod_positions.append(rod_positions[-1]+diameter_rod/2.+shift)
        center = (rod_positions[0]+rod_positions[-1])/2.
        
        rigid.set_param('rod',
                        types=['r']+['R']*rod_segments+['s'],
                        diameters=[2*radius_sticker]+[diameter_rod]*rod_segments+[2*radius_sticker],
                        positions=[(x-center,0,0) for x in rod_positions],
                        )
        
        linker_shift = radius_sticker+diameter_linker/2.0
        rigid.set_param('linker',
                        types=['l']+['s']*5,
                        diameters=[2*radius_sticker]*6,
                        positions=[ (-linker_shift,0,0), (linker_shift,0,0),
                                    (0,-linker_shift,0), (0,linker_shift,0),
                                    (0,0,-linker_shift), (0,0,linker_shift),
                                  ],
                        )
        rigid.create_bodies()
        snapshot = system.take_snapshot()
    else:
        #don't set up system, load from file to continue
        if not os.path.exists(inputfile):
            print("Error: inputfile %s does not exist"%inputfile)
            sys.exit(1)
        system=hoomd.init.read_gsd(inputfile,frame=-1)
        rigid = hoomd.md.constrain.rigid();
        rod_positions = [0]
        for i in range(rod_segments+1):
            if i == 0 or i == rod_segments: 
                shift = radius_sticker
            else:
                shift = 0
            rod_positions.append(rod_positions[-1]+diameter_rod/2.+shift)
        center = (rod_positions[0]+rod_positions[-1])/2.
        
        rigid.set_param('rod',
                        types=['r']+['R']*rod_segments+['s'],
                        diameters=[2*radius_sticker]+[diameter_rod]*rod_segments+[2*radius_sticker],
                        positions=[(x-center,0,0) for x in rod_positions],
                        )
        
        linker_shift = radius_sticker+diameter_linker/2.0
        rigid.set_param('linker',
                        types=['l']+['s']*5,
                        diameters=[2*radius_sticker]*6,
                        positions=[ (-linker_shift,0,0), (linker_shift,0,0),
                                    (0,-linker_shift,0), (0,linker_shift,0),
                                    (0,0,-linker_shift), (0,0,linker_shift),
                                  ],
                        )
        rigid.validate_bodies()
        snapshot = system.take_snapshot()
        crowder_type = snapshot.particles.types.index('crowder')
        number_crowders = (snapshot.particles.typeid == crowder_type).sum()
        
    
    nl = md.nlist.cell()
    nl.reset_exclusions(exclusions = ['1-2', 'body','constraint'])
    
    #pair = md.pair.gauss(r_cut=diameter_crowder, nlist=nl)
    #pair.set_params(mode='shift')
    #alpha=0 is pure repulsive
    #pair = md.pair.force_shifted_lj(r_cut=diameter_crowder, nlist=nl)
    pair = md.pair.table(width=1000,nlist=nl)
    radius_dict = {
#should this be bigger to prevent overlaps?
                    #'R': diameter_rod/2.,
                    'R': diameter_rod/2*1.2,
                    'r': radius_sticker,
                    'linker': diameter_linker/2.,
                    'rod': diameter_rod/2.,
                    'crowder': diameter_crowder/2.,
                    'l': radius_sticker,
                    's': radius_sticker,
    }

    if lj_epsilon > 0 :
        pair_lj = md.pair.lj(r_cut=diameter_crowder,nlist=nl)
        pair_lj.set_params(mode="xplor")
    
    for p_type_1 in radius_dict.keys():
        for p_type_2 in radius_dict.keys():
            sigma = radius_dict[p_type_1] + radius_dict[p_type_2]
            #if p_type_1 == 'rod' or p_type_2 =='rod' or p_type_1 == 'l' or p_type_2 == 'l':
            if p_type_1 == 'R' and p_type_2 == 'R':
                epsilon = rod_repulsion
            else:
                epsilon = sphere_repulsion

            if p_type_1 in ['l','r','s'] or p_type_2 in ['l','r','s']:
                epsilon = 0

            
            pair.pair_coeff.set(p_type_1, p_type_2, func=softrepulsion, rmin=0.*sigma,rmax=1.5*sigma,coeff=dict(epsilon=epsilon,rcut=sigma))
            if lj_epsilon>0:
                if p_type_1 in ["linker","rod"] and p_type_2 in ["linker","rod"]:
                    pair_lj.pair_coeff.set(p_type_1,p_type_2,epsilon=lj_epsilon,sigma=sigma,r_cut = 2.5*sigma)
                else:
                    pair_lj.pair_coeff.set(p_type_1,p_type_2,epsilon=0,sigma=sigma,r_cut = -2.5*sigma)
                
            
    ##DYNAMIC BONDING PART:
    #
    particle_clusterids = []
    cluster_id = 0
    for i in range(len(snapshot.particles.typeid)):
        particle_clusterids.append("%i"%i)
    particle_clusterids = " ".join(particle_clusterids)
    
    checksteps=10
    kT=1.0
    r0=bond_length
    #standarddev=np.sqrt(kT/spring_constant)*r0
    standarddev=np.sqrt(kT/spring_constant)

    binding_distance_in_std = binding_distance/standarddev
    
    harmonic=hoomd.md.bond.harmonic()
    harmonic.bond_coeff.set('r-l',k=spring_constant,r0=r0)
    rmax=r0+standarddev*binding_distance_in_std
    rmin=r0-standarddev*binding_distance_in_std
    print("Bonding if between %f and %f"%(rmin,rmax))
    #
    Tmelt=50.0
    #alpha=0.05
    alpha=10
    flag_temperature=0
    flag_metropolis=1
    flag_force=0
    #
    kon_init=1/(dt*checksteps)
    kon_melt=0
    #
    koff_init=koff         #koff_init set to 0 

    koff_melt = kon_init + kon_melt - koff_init
    #
    permanent_bonds = 0
    
    group_r = hoomd.group.type(name='r-particles', type='r')
    group_l = hoomd.group.type(name='l-particles', type='l')
    group_bondable = hoomd.group.union('bondable',group_r,group_l)
    
    
    updater = db.update.dybond(nl, group=hoomd.group.all(), period=checksteps)
    updater.set_params(bond_type='r-l',A='r',B='l',nondybonds=permanent_bonds,r0=r0,rmin=rmin,rmax=rmax,kspring=spring_constant,flag_temperature=flag_temperature,flag_force=flag_force,flag_metropolis=flag_metropolis,Tmelt=Tmelt,alpha=alpha,kon_init=kon_init,kon_melt=kon_melt,koff_init=koff_init,koff_melt=koff_melt,checksteps=checksteps,dt=dt,particle_clusterids=particle_clusterids,self_avoiding=0,userseed=seed) 

#    #first minimize
#    fire=md.integrate.mode_minimize_fire(dt=dt/10, ftol=1e-2, Etol=1e-7)
#    nve=md.integrate.nve(group=hoomd.group.rigid_center())
#    while not(fire.has_converged()):
#       hoomd.run(5000)
#    nve.disable()
    
    ##integrate at constant temperature
    md.integrate.mode_standard(dt=dt)
    hoomd.md.integrate.langevin(group=hoomd.group.rigid_center(), kT=kT, seed=seed, dscale=0.001)
    if number_crowders > 0:
        hoomd.md.integrate.langevin(group=hoomd.group.nonrigid(), kT=crowder_temperature, seed=seed, dscale=0.001)
    #
    hoomd.dump.gsd(gsd_file,period=dump_frames,group=hoomd.group.all(),overwrite=True,
           #dynamic=['attribute','momentum','topology'])
           dynamic=['topology'])
    #
    logger=hoomd.analyze.log(log_file,quantities=['potential_energy', 'bond_harmonic_energy','temperature'],period=dump_frames,overwrite=True)
    #
    
    hoomd.run(nsteps)
