import hoomd
import numpy as np
import os
import sys
import argparse
import hoomd.dybond_plugin as db
from hoomd import md
from potentials import softrepulsion
from run_packmol import gen_lattice_packmol
#import numba as nb

#@nb.jit(nopython=True)
def wca(epsilon, sigmas, distances):
    wca_cuts = 2**(1./6)*sigmas
    sr = sigmas/distances
    energy = 4*epsilon*((sr)**12-(sr)**6) + epsilon
    energy[distances>wca_cuts] = 0
    return np.sum(energy)

#@nb.jit(nopython=True)
def gen_lattice_mc(number_particles,box_length,particle_type_list, diameter_list, max_iter=20000,n_attempts=3, temperature=1.1, harmonic_const=10):
    success = False
    min_position = -box_length/2.
    box_size = np.array([box_length,box_length,box_length])

    for att in range(n_attempts):
        particle_position_array = np.zeros((number_particles,3))
        radius_list = diameter_list[particle_type_list]/2.
        prev_energy = 0
        for i in range(number_particles-1,-1,-1):
            radius_i = diameter_list[particle_type_list[i]]/2.
            overlap_dists = radius_list+radius_i
            for j in range(max_iter):
                particle_position_attempt = min_position*np.ones(3) + np.random.random(size=3)*box_length
                if i == number_particles-1:
                    particle_position_array[i,:] = particle_position_attempt
                    break
                else:
                    dr = particle_position_array - particle_position_attempt
                    dr = dr - box_size*np.floor(dr/box_size+0.5)
                    dists = np.sqrt( (dr*dr).sum(axis=1) )

                    if temperature == 0:
                        if (dists<overlap_dists[:len(dists)]).sum()>0:
                            continue
                        else:
                            particle_position_array[i,:] = particle_position_attempt
                            break
                    else:
                        energy = harmonic_const*np.sum( (dists-overlap_dists[:len(dists)])**2 )
                        #energy = wca(harmonic_const, overlap_dists[:len(dists)], dists)
                        dE = energy-prev_energy
                        if np.random.random()<np.exp(dE/temperature):
                            particle_position_array[i,:] = particle_position_attempt
                            prev_energy = energy
                            break
                        else:
                            continue

            if j >= max_iter -1:
                print("Failed to add particle on", max_iter," attemps for particle ",i,"on attempt",att)
                break
        if i==0:
            success = True
            break
    if att == n_attempts-1 and success is False:
        print("No mc solution found")
        return None

    return particle_position_array

def gen_lattice_v1(number_density,num_atoms,fudge=0.95):
    box_length = (number_density/num_atoms)**(-1./3)
    repeats = int(np.ceil((num_atoms/2.)**(1./3)))
    num_positions = 2*(repeats**3)

    lattice_spacing = 1.22*(number_density/num_positions*num_atoms)**(-1./3)

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
    positions = positions/(positions.max(axis=0)-positions.min(axis=0))*box_length*fudge
    positions -= box_length/2. * np.ones(positions.shape) #- lattice_spacing/2.0*np.ones(positions.shape)
    return positions


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--gpu",default=False,action="store_true")
    parser.add_argument("--outprefix",help="Start of output simulation name (default: %(default)s)",default=None,type=str,required=True)
    parser.add_argument("--nsteps",help="Number of steps of simulation to run (default: %(default)s)",default=None,type=int,required=True)
    parser.add_argument("--dt",help="Simulation time step (default: %(default)s)",default=0.002,type=float)
    parser.add_argument("--koff",help="Off rate for unbinding in time units (default: %(default)s)",default=0.001,type=float)
    parser.add_argument("--sphere_repulsion",help="Repulsion of spheres including crowders in units of kT, for soft potential (default: %(default)s)",default=500,type=float)
    parser.add_argument("--rod_repulsion",help="Repulsion of spheres including crowders in units of kT, for soft potential (default: %(default)s)",default=500,type=float)
#box_length=860.
    parser.add_argument("--crowder_temperature",help="Temperature of crowders, relative to 1.0 (default: %(default)s)",default=1.0,type=float)
    parser.add_argument('--initial_box_length',help='Initial side length of simulation box in nm, for compression',default=None,type=float)
    parser.add_argument('--compression_steps',help='Number of steps to run compression for',default=500000,type=int)
    parser.add_argument("--box_length",help="Side length of simulation box in nm (default: %(default)s)",default=860,type=float)
    parser.add_argument("--binding_distance",help="Maximum distance where binding can take place, in nm (default: %(default)s)",default=0.9,type=float)
    parser.add_argument("--number_rods",help="Number of rod proteins (default: %(default)s)",default=1170,type=int)
    parser.add_argument("--number_linkers",help="Number of rod proteins (default: %(default)s)",default=390,type=int)
    parser.add_argument("--number_gems",help="Number of GEM proteins (default: %(default)s)",default=0,type=int)
    parser.add_argument("--volume_fraction_ribosome",help="Fraction of volume taken up by 30nm ribosome (default: %(default)s)",default=0,type=float)
    parser.add_argument("--volume_fraction_polysome",help="Fraction of volume taken up by 100nm polysome (default: %(default)s)",default=0,type=float)
    parser.add_argument("--seed",help="Random seed (default: %(default)s)",default=1,type=int)
    parser.add_argument("--dump_frequency",help="Dump frequency in time units(default: %(default)s)",type=int,default=50)
    parser.add_argument("--lj_epsilon",help="Epsilon for LJ attraction (default: %(default))",type=float,default=0)
    parser.add_argument("-i","--inputfile",default=None,help="Input file for continuing (gsd file, optional)",type=str)
    parser.add_argument('--closed_box',default=False,action="store_true",help="Turn on enclosing walls")
#    parser.add_argument('--write_unwrapped',default=False,action="store_true",help="Write unwrapped trajectory")
    parser.add_argument("--wall_epsilon",help="Epsilon for wall repulsion (default: %(default))",type=float,default=0.1)
    parser.add_argument("--wall_shift",help="Wall shift in nm to prevent explosion (default: %(default))",type=float,default=10)
    parser.add_argument('--diameter_gem',help="Diameter of gem particle (default: %(default))",type=float,default=40)
    parser.add_argument('--gamma_scale',help="Friction coefficient to multiply diameter to give diffusion coeff (default: %(default))",type=float,default=0.001)

    args = parser.parse_args()
    print("Command line arguments given:")
    print(args)

    locals().update(vars(args))
    if gpu:
        hoomd.context.initialize("--mode=gpu")
        print("Running on the GPU")
    else:
        hoomd.context.initialize("--mode=cpu")

    print("Box length: ",box_length)

    gsd_file = outprefix+".gsd"
    gsd_file_compress = outprefix+".compress.gsd"

    gsd_file_unwrapped = outprefix+".unwrapped.gsd"
    log_file = outprefix+".log"
    out_dir = os.path.dirname(gsd_file)

    np.random.seed(seed)
    dump_frames = int(dump_frequency/dt)

    rod_segments = 2
    diameter_rod = 23.4/float(rod_segments)
    diameter_linker = 12.6
    diameter_ribosome = 30.
    diameter_polysome = 100.
    diameter_gem = args.diameter_gem
    radius_sticker = 1
    bond_length=2*radius_sticker
    spring_constant=20

    radius_dict = {
#should this be bigger to prevent overlaps?
                    'R': (diameter_rod/2)*1.15,   
 #                   'R': diameter_rod/2,
                    'r': radius_sticker,
                    'linker': diameter_linker/2.,
                    'rod': diameter_rod/2.,
                    'ribosome': diameter_ribosome/2.,
                    'polysome': diameter_polysome/2.,
                    'gem': diameter_gem/2.,
                    'l': radius_sticker,
                    's': radius_sticker,
    }

    print("Radii of the particles:")
    print(radius_dict)

    if inputfile is None:
        volume_in_nm = (box_length)**3
        volume_in_microns = (box_length/1000.)**3
        if initial_box_length is not None:
            volume_in_microns_initial = (initial_box_length/1000.)**3
    
    #particles / micron^3
    #number_density_linkers=1839
    #number_density_rods=5517
    #number_density_linkers = 609 
    #number_density_rods = 1827
    #number_linkers = int(volume_in_microns*number_density_linkers)
    #number_rods = int(volume_in_microns*number_density_rods)
    
        #in nm 3
        volume_ribosome = 4./3*np.pi*((diameter_ribosome/2.)**3)
        volume_polysome = 4./3*np.pi*((diameter_polysome/2.)**3)
        number_ribosomes = int(volume_fraction_ribosome*volume_in_nm/volume_ribosome)
        number_polysomes = int(volume_fraction_polysome*volume_in_nm/volume_polysome)
        num_per_type = [ number_rods, number_linkers , number_gems , number_ribosomes , number_polysomes]
        particle_type_list = np.repeat(np.arange(len(num_per_type)).astype(int),repeats=num_per_type)
        number_particles = np.sum(num_per_type)
        print("Starting the system with %i rods , %i linkers, %i gems,  %i ribosomes, %i polysomes"%tuple(num_per_type))
    
        particle_types = particle_types=['rod', 'linker','gem', 'ribosome','polysome','r','l','s','R']
        #diameter_list = np.array([diameter_rod, diameter_linker, diameter_gem, diameter_ribosome, diameter_polysome ])
        diameter_list = [radius_dict[i]*2 for i in particle_types]
        diameter_list_pack = np.array([diameter_rod*3, diameter_linker*1.1, diameter_gem, diameter_ribosome, diameter_polysome ])

        scale_factor = 1.0
        if closed_box is True:
            if initial_box_length is not None:
                print("Warning: changing box size may not yet be supported with walls!")
            full_box_length = box_length + 4*diameter_list.max()
            print("Using a closed box")
        elif initial_box_length is not None:
            full_box_length = initial_box_length
            density_scale_factor = (initial_box_length/box_length)**-3
        else:
            full_box_length = box_length

        snapshot = hoomd.data.make_snapshot(N=number_particles, box=hoomd.data.boxdim(L=full_box_length,dimensions=3), particle_types=particle_types, 
                      bond_types=['r-l'])


        number_density = number_particles / volume_in_nm
        print("Generating initial configuration")
        random_lattice_positions = gen_lattice_v1(number_density*density_scale_factor,number_particles)
        #random_lattice_positions = gen_lattice_mc(number_particles,box_length,particle_type_list, diameter_list, harmonic_const=sphere_repulsion)#, seed=seed)
        #random_lattice_positions = gen_lattice_packmol(box_length,particle_types, particle_type_list, diameter_list_pack, outprefix=outprefix)
        if random_lattice_positions is None:
            print("Error, failed to generate starting confuration by MC")
            sys.exit(1)
        print("min and max positions:")
        print(random_lattice_positions.min(axis=0))
        print(random_lattice_positions.max(axis=0))
        
        for i in range(number_particles):
            typeid = particle_type_list[i]
            snapshot.particles.typeid[i] = typeid
            #note, diameter of R is set in the rigid code above
            snapshot.particles.diameter[i]=diameter_list[typeid]
            if typeid == 0: #rod
                snapshot.particles.moment_inertia[i] = np.array([1,1,0])*0.1
            else:
                snapshot.particles.moment_inertia[i] = np.array([1,1,1])*0.1
            snapshot.particles.mass[i]=1
            snapshot.particles.position[i] = random_lattice_positions[i]
        
        snapshot.bonds.resize(0)
        #snapshot.bonds.group[0]=[0,1]
        system=hoomd.init.read_snapshot(snapshot)
        
#        # Add constituent particles of type rod and create the rods
#        system.particles.types.add('R');
        
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
                        diameters=[2*radius_sticker]+[2*radius_dict['R']]*rod_segments+[2*radius_sticker],
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
        ribosome_type = snapshot.particles.types.index('ribosome')
        number_ribosomes = (snapshot.particles.typeid == ribosome_type).sum()
        polysome_type = snapshot.particles.types.index('polysome')
        number_polysomes = (snapshot.particles.typeid == polysome_type).sum()
        gem_type = snapshot.particles.types.index('gem')
        number_gems = (snapshot.particles.typeid == gem_type).sum()
        
    
    nl = md.nlist.cell()
    nl.reset_exclusions(exclusions = ['body','constraint'])
    
    #pair = md.pair.gauss(r_cut=diameter_crowder, nlist=nl)
    #pair.set_params(mode='shift')
    #alpha=0 is pure repulsive
    #pair = md.pair.force_shifted_lj(r_cut=diameter_crowder, nlist=nl)

    #use this for soft repulsion
    pair = md.pair.table(width=1000,nlist=nl)
    #pair_wca = md.pair.lj(r_cut=diameter_ribosome,nlist=nl)
    #pair_wca.set_params(mode="shift")


    if lj_epsilon > 0 :
        pair_lj = md.pair.lj(r_cut=diameter_ribosome,nlist=nl)
        pair_lj.set_params(mode="xplor")

    #if closed, turn on walls
    if closed_box is True:
        wall_locations = np.array([ [-box_length/2-wall_shift,0,0], \
                           [+box_length/2+wall_shift,0,0], \
                           [0,-box_length/2-wall_shift,0], \
                           [0,+box_length/2+wall_shift,0], \
                           [0,0,-box_length/2-wall_shift], \
                           [0,0,+box_length/2+wall_shift], \
                        ])
        wall_group = hoomd.md.wall.group()
        for l in wall_locations:
            wall_group.add_plane(l,-l/np.linalg.norm(l), inside=True)
        wall_repulsion=hoomd.md.wall.lj(wall_group, r_cut=2.5*np.max(diameter_list))
    
    for p_type_1 in radius_dict.keys():
        if closed_box is True:
            wall_repulsion.force_coeff.set(p_type_1,sigma=radius_dict[p_type_1]/3.,epsilon=wall_epsilon,alpha=0,r_cut=2.5*radius_dict[p_type_1])

        for p_type_2 in radius_dict.keys():
            sigma = radius_dict[p_type_1] + radius_dict[p_type_2]
            #if p_type_1 == 'rod' or p_type_2 =='rod' or p_type_1 == 'l' or p_type_2 == 'l':
            if p_type_1 in ('R','rod','linker') and p_type_2 in ('R','rod','linker'):
                epsilon = rod_repulsion
            else:
                epsilon = sphere_repulsion

            if p_type_1 in ['l','r','s'] and p_type_2 in ['l','r','s']:
                epsilon = 0
            
            pair.pair_coeff.set(p_type_1, p_type_2, func=softrepulsion, rmin=0.*sigma,rmax=1.5*sigma,coeff=dict(epsilon=epsilon,rcut=sigma))
            #pair_wca.pair_coeff.set(p_type_1, p_type_2,epsilon=epsilon,sigma=sigma,r_cut=(2**(1./6)*sigma))
            #pair_wca.pair_coeff.set(p_type_1, p_type_2,epsilon=epsilon,sigma=sigma,r_cut=sigma)

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
        #old, everything unique
        #particle_clusterids.append("%i"%i)
        #new, by body
        body_code = snapshot.particles.body[i]
        if body_code > 1e5: 
            body_code=-1
        particle_clusterids.append("%i"%body_code)
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
    flag_force=0
    flag_metropolis=1
    #
    kon_init=1/(dt*checksteps)
    kon_melt=0
    #
    koff_init=koff         #koff_init set to 0 

    koff_melt = kon_init + kon_melt - koff_init
    #
    #permanent_bonds = snapshot.bonds.N
    permanent_bonds = 0
    
    group_r = hoomd.group.type(name='r-particles', type='r')
    group_l = hoomd.group.type(name='l-particles', type='l')
    group_bondable = hoomd.group.union('bondable',group_r,group_l)
    
    
    
    group_gem = hoomd.group.type(name='gem-particles', type='gem')
    brownian_group = hoomd.group.union(name='move',a=hoomd.group.rigid_center(),b=group_gem)
    if number_ribosomes > 0 or number_gems>0 or number_polysomes>0:
        group_polysome = hoomd.group.type(name='polysome-particles', type='polysome')
        group_ribosome = hoomd.group.type(name='ribosome-particles', type='ribosome')
        group_crowder = hoomd.group.union('crowder',group_ribosome,group_polysome)
        min_group = hoomd.group.union(name='min',a=brownian_group,b=group_crowder)
    else:
        min_group = brownian_group

    if inputfile is None:
        #first minimize
        fire=md.integrate.mode_minimize_fire(dt=dt*0.3, ftol=1e-3, Etol=1e-7)
        nve=md.integrate.nve(group=min_group)
        while not(fire.has_converged()):
           hoomd.run(20000)
        nve.disable()

  ##integrate at constant temperature
    md.integrate.mode_standard(dt=dt)
    #ld1 = hoomd.md.integrate.brownian(group=brownian_group, kT=kT, seed=seed, dscale=gamma_scale)
    ld1 = hoomd.md.integrate.langevin(group=brownian_group, kT=kT, seed=seed, dscale=gamma_scale)
    if number_ribosomes > 0 or number_gems>0 or number_polysomes>0:
        if number_ribosomes+number_polysomes>0:
            #ld2 = hoomd.md.integrate.brownian(group=group_crowder, kT=crowder_temperature, seed=seed, dscale=gamma_scale)
            ld2 = hoomd.md.integrate.langevin(group=group_crowder, kT=crowder_temperature, seed=seed, dscale=gamma_scale)

    if initial_box_length is not None and inputfile is None:
        print(f"Compressing from box length {initial_box_length} to {box_length} in {compression_steps} steps")
        compress_dump = hoomd.dump.gsd(gsd_file_compress,period=dump_frames//5,group=hoomd.group.all(),overwrite=True,dynamic=['topology','image'])
        compress_logger=hoomd.analyze.log(log_file.replace('.log','.compress.log'),quantities=['potential_energy', 'bond_harmonic_energy','temperature'],period=dump_frames,overwrite=True)

        # compress the box before starting and enabling bonding
        L = hoomd.variant.linear_interp(points = [(0,initial_box_length), (compression_steps,box_length)],zero=0)
        hoomd.update.box_resize(L=L,period=1,phase=0)
        hoomd.run(compression_steps)

        compress_dump.disable()
        compress_logger.disable()

    updater = db.update.dybond(nl, group=hoomd.group.all(), period=checksteps)
    updater.set_params(bond_type='r-l',A='r',B='l',nondybonds=permanent_bonds,r0=r0,rmin=rmin,rmax=rmax,kspring=spring_constant,flag_temperature=flag_temperature,flag_force=flag_force,flag_metropolis=flag_metropolis,Tmelt=Tmelt,alpha=alpha,kon_init=kon_init,kon_melt=kon_melt,koff_init=koff_init,koff_melt=koff_melt,checksteps=checksteps,dt=dt,particle_clusterids=particle_clusterids,self_avoiding=0,userseed=seed) 
     
    hoomd.dump.gsd(gsd_file,period=dump_frames,group=hoomd.group.all(),overwrite=True,
           #dynamic=['attribute','momentum','topology'])
           dynamic=['topology','image'])
    
    if closed_box is True:
        logger=hoomd.analyze.log(log_file,quantities=['potential_energy', 'bond_harmonic_energy','external_wall_lj_energy','temperature'],period=dump_frames,overwrite=True)
    else:
        logger=hoomd.analyze.log(log_file,quantities=['potential_energy', 'bond_harmonic_energy','temperature'],period=dump_frames,overwrite=True)
        
    hoomd.run(nsteps)
