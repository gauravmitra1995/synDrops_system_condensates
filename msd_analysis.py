
#simple version
def compute_msd_xyz_data(xyz):
    import numpy as np
    msd_avg = []
    msd_std = []
    nframes = len(xyz)

    for frame_sep in range(1,nframes-1):
        tmp_msd = []
        for start_frame in range(0,nframes-frame_sep):
            dr = xyz[start_frame+frame_sep]-xyz[start_frame]
            dr2 = np.sum(dr*dr,axis=-1)
            tmp_msd.extend(dr2)
        msd_avg.append(np.mean(tmp_msd))
        msd_std.append(np.std(tmp_msd))
 
    return np.array(msd_avg),np.array(msd_std)

def msd_trajectory(trajectory_file,dt,typeid, first_frame=0):
    from sklearn.linear_model import LinearRegression
    import numpy as np
    from gsd import hoomd as gsd
    import matplotlib.pyplot as plt
    
    trajectory = gsd.open(trajectory_file,'rb') # read gsd file
    coordinates_type0 = []
    
    frame_stepnums = []
    prev_cfg = None
    total_shift = None

    print("# Analyzing for type:",trajectory[0].particles.types[typeid])

    frame_list = [frame for frameidx,frame in enumerate(trajectory) if frameidx>first_frame]

    for frame in frame_list:
        cfg = frame.configuration
        frame_stepnums.append(cfg.step)
        particles = frame.particles
        particle_types = particles.typeid
        center_positions = particles.position[particle_types==typeid,:]
        if total_shift is None :
            total_shift = np.zeros(np.shape(center_positions))
        
        #xyz length of simulation box:
        box = cfg.box[:cfg.dimensions]
        #compute shift from crossing periodic box
        #shift = particles.image[particle_types==1]*box

        if prev_cfg is not None:
            total_shift -= box*np.floor( (center_positions - prev_cfg)/box + 0.5 )
        prev_cfg = np.copy(center_positions)
        center_positions = center_positions+total_shift
        coordinates_type0.append(center_positions)
    #to debug wall, make sure there is no shift
    #print(total_shift)
    
    coordinates_type0 = np.array(coordinates_type0)
#    print(coordinates_type0[:,1,0])
#    plt.plot(coordinates_type0[:,:,0],c='r')
#    plt.plot(coordinates_type0[:,:,1],c='b')
#    plt.plot(coordinates_type0[:,:,2],c='k')
#    plt.show()
#    sys.exit()
    msd, msd_std = compute_msd_xyz_data(coordinates_type0)
    #plt.errorbar(frame_stepnums[:-2],msd,msd_std,linestyle='-')
    times = dt*np.array(frame_stepnums[:-2])
    plt.plot(times,msd,linestyle='-')
    plt.ylabel('MSD(t)')
    plt.xlabel('time (HOOMD units)')
    #plt.yscale('log')
    #plt.xscale('log')

    lr = LinearRegression()
#    lr.fit(times.reshape(-1,1),msd.reshape(-1,1))
#    msd_predict = lr.predict(times.reshape(-1,1))
#    D = float(lr.coef_/6)

    #log fit
    #fit_slice_length = len(times)//10
    fit_slice_length = len(times)
    lr.fit(np.log10(times.reshape(-1,1))[:fit_slice_length],np.log10(msd.reshape(-1,1)[:fit_slice_length]))
    msd_predict = lr.predict(np.log10(times.reshape(-1,1)))
    D=(10**lr.intercept_)/6
    alpha = lr.coef_[0]

    print("Fit D=%f nm^2/t"%(D))
    print("Fit alpha=%f"%(alpha))
    print("Time unit t = %e s"%(D/3e5)) # ~1micron^2/s = 1e6 nm^2/s, 0.3 micron^2/s=3e5nm^2/s
    plt.plot(times,10**msd_predict,linestyle='--',label="D=%.3f,alpha=%.3f"%(D,alpha))
    plt.yscale('log')
    plt.xscale('log')

    plt.legend()
    plt.show()

if __name__ == "__main__":
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument("--file",type=str,required=True)
    parser.add_argument("--typeid",type=int,default=2)
    parser.add_argument("--dt",type=float,default=0.002)
    parser.add_argument("--first_frame",type=int,default=0)
    args = parser.parse_args()
    msd_trajectory(args.file,args.dt,args.typeid,args.first_frame)
