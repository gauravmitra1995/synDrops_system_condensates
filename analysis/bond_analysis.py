import sys
import numpy as np 

def gem_analysis(trajectory,frame_id=-1,step_time=1,gem_type='gem'):
    frame = trajectory[frame_id]
    box_size = np.array(frame.configuration.box.flat)[:3]
    step = frame.configuration.step
    bodies = frame.particles.body
    types = frame.particles.typeid

    type_labels = frame.particles.types
    gems = np.where( types==type_labels.index(gem_type) )[0]

    all_xyz = frame.particles.position
    main_particle_xyz = all_xyz[gems]

    if step > 0:
        frame_prev = trajectory[frame_id-1]
        step_prev = frame_prev.configuration.step
        delta_t = (step - step_prev)*step_time
        main_particle_xyz_prev = frame_prev.particles.position[gems]

        dx = main_particle_xyz - main_particle_xyz_prev
        dx = dx - box_size * np.floor(dx/box_size + 0.5)
        dx2 = (dx*dx).sum(axis=-1)
        diff = dx2/(6*delta_t)/1e6 #divide by a million to convert to microns squared
    else:
        diff = np.zeros( len(main_particle_xyz))
        dx = None

    position_array = np.append(main_particle_xyz, diff.reshape(-1,1),axis=-1)
    indices = np.arange(len(main_particle_xyz)).reshape(-1,1)

    position_array = np.append(indices,position_array,axis=-1)

    return gems, position_array, dx

def get_bond_table(trajectory,frame_id=-1,sort_type="hexamer", step_time=1):
    if sort_type == "hexamer":
        sort_column = 1
    else:
        sort_column = 0

    frame = trajectory[frame_id]

    bodies = frame.particles.body

    bonds = frame.bonds
    bond_tags = bonds.group

    bond_tags = bond_tags[np.argsort(bond_tags[:,sort_column])]
    
    column1_bodies = bodies[bond_tags[:,0]]
    column2_bodies = bodies[bond_tags[:,1]]

    body_array = np.append(column1_bodies.reshape(-1,1),column2_bodies.reshape(-1,1),axis=-1)

    return body_array,bodies

def bonds_analysis(trajectory,frame_id=-1,sort_type="hexamer", step_time=1, notebook=True):
    if sort_type == "hexamer":
        sort_column = 1
    else:
        sort_column = 0
    frame = trajectory[frame_id]
    box_size = np.array(frame.configuration.box.flat)[:3]
    step = frame.configuration.step
    bodies = frame.particles.body
    types = frame.particles.typeid
    bonds = frame.bonds
    bond_tags = bonds.group
    bond_tags = bond_tags[np.argsort(bond_tags[:,sort_column])]

    #is this needed? as long as bonds end up being unique?
    #if len(bond_tags)>0:
    #    assert np.all(types[bond_tags[:,0]]==types[bond_tags[0,0]]),"Not all particles in list are the same type"
    #    assert np.all(types[bond_tags[:,1]]==types[bond_tags[0,1]]),"Not all particles in list are the same type"
    column1_bodies = bodies[bond_tags[:,0]]
    column2_bodies = bodies[bond_tags[:,1]]

    body_array = np.append(column1_bodies.reshape(-1,1),column2_bodies.reshape(-1,1),axis=-1)

    all_xyz = frame.particles.position
    main_particle_xyz = all_xyz[types < 3]
    #original_indices = np.arange(len(all_xyz)).astype(int)[types<3].reshape(-1,1)

    if step > 0:
        frame_prev = trajectory[frame_id-1]
        step_prev = frame_prev.configuration.step
        delta_t = (step - step_prev)*step_time
        main_particle_xyz_prev = frame_prev.particles.position[types<3]

        dx = main_particle_xyz - main_particle_xyz_prev
        dx = dx - box_size * np.floor(dx/box_size + 0.5)
        dx2 = (dx*dx).sum(axis=-1)
        diff = dx2/(6*delta_t)/1e6 #divide by a million to convert to microns squared
    else:
        diff = np.zeros( len(main_particle_xyz))

    position_array = np.append(main_particle_xyz, diff.reshape(-1,1),axis=-1)
    indices = np.arange(len(main_particle_xyz)).reshape(-1,1)

    if notebook is True:
        return body_array, position_array, diff
    else:
        position_array = np.append(indices,position_array,axis=-1)
        #position_array = np.append(position_array,original_indices,axis=-1)
        return position_array

def nbonds_v_time(trajectory,bonding_type=5):
    max_bonds = (trajectory[0].particles.typeid == bonding_type).sum() 
    n_bonds = []
    for i in range(len(trajectory)):
        n_bonds.append(len(bonds_analysis(trajectory,i)[0]))
    n_bond_array = np.array(n_bonds)
    frac_bond_array = np.array(n_bonds).astype(float)/max_bonds
    result_array = np.append(n_bond_array.reshape(-1,1),frac_bond_array.reshape(-1,1),axis=-1)
    return result_array
    
if __name__ == "__main__":
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument("file",type=str)
    parser.add_argument("--frame",type=int,required=False,default=-1)
    parser.add_argument("--trajectory",default=False,action="store_true")
    parser.add_argument("--binding_const",default=False,action="store_true")
    parser.add_argument("--step_time",default=0.134,type=float,help="Time of an MD step in microseconds (default %(default)d)")
    parser.add_argument("--analysis_stride",default=10,type=int,help="How many frames to skip in writing bond analysis to a file (default %(default)d)")
    parser.add_argument("--output_prefix",default=None)
    args = parser.parse_args()

    from gsd import hoomd as gsd
    trajectory = gsd.open(args.file,'rb') # read gsd file
    if args.trajectory is True:
        nbonds_trajectory = nbonds_v_time(trajectory)
        np.savetxt(sys.stdout,nbonds_trajectory.reshape(-1,2),fmt="%d %.2f") 
    elif args.binding_const is True:
        box_size = np.array(trajectory[0].configuration.box.flat)[:3]
        volume_in_nm3 = np.prod(box_size)
        #just use last 1% of trajectory
        nbonds_trajectory = nbonds_v_time(trajectory)[len(trajectory)//99:]
        num_particles = len(np.unique(trajectory[0].particles.body))
        num_a_total = num_particles/2.
        num_ab = np.mean(nbonds_trajectory[:,0])
        frac_bonded = np.mean(nbonds_trajectory[:,1])

        kd = (num_a_total - num_ab)**2/(num_ab * volume_in_nm3) * 10/6.

        print("Frac bonded:",frac_bonded)
        print("Kd(M):",kd)
        print("Kd(muM):",kd*1e6)
    else:
        #put step time in to seconds
        if args.output_prefix is None:
            bonded_particles, particle_positions = bonds_analysis(trajectory,frame_id=args.frame,step_time=args.step_time/1e6)
            np.savetxt(sys.stdout,bonded_particles,fmt='%d',delimiter=',')
        else:
            for i in range(0,len(trajectory),args.analysis_stride):
                step = trajectory[i].configuration.step
                step_time_seconds = step*args.step_time/1e6
                bonded_particles, particle_positions = bonds_analysis(trajectory,frame_id=i,step_time=args.step_time/1e6)

                np.savetxt(args.output_prefix+'_%.3f.nLLinks'%step_time_seconds,bonded_particles,fmt='%d',delimiter=',')
                np.savetxt(args.output_prefix+'_%.3f.nLCoords'%step_time_seconds,particle_positions,fmt='%d,%.3e,%.3e,%.3e,%.3e',delimiter=',')
