#source https://04306371071604527041.googlegroups.com/attach/1b66df8c50f03/msd.py?part=0.1.1&view=1&vt=ANaJVrHUo20ajDc7GziSA4NQDVbGy2ZwerDjw7XmY8_YMI93742W1a0t0I2CrgwMJxzIEW-eD_YrwyAbqJ8Tk-xdimG_BLSKC5basnvwGuJgeXa9CYI4iJo
from gsd import fl
from gsd import hoomd as gsdhoomd
from freud.box import Box
import numpy as np
import sys

def autocorrFFT_n(x):
  N=len(x)
  F = np.fft.fft(x, n=2*N, axis=0)  #2*N because of zero-padding
  PSD = F * F.conjugate()
  res = np.fft.ifft(PSD, axis=0)
  res= (res[:N]).real   #now we have the autocorrelation in convention B
  n=N*np.ones(N)-np.arange(0,N) #divide res(m) by (N-m)
  return res/n[:,None] #this is the autocorrelation in convention A

def msd_fft_n(r):
  if r.ndim == 2:
      r = np.asarray(r).reshape((-1,1,3))
  N=len(r)
  D=np.square(r).sum(axis=2)
  D=np.append(D,np.zeros(r.shape[:2]),axis=0)
  S2=np.sum([autocorrFFT_n(r[:,:, i]) for i in range(r.shape[2])],axis=0)
  Q=2*D.sum(axis=0)
  S1=np.zeros(r.shape[:2])
  for m in range(N):
      Q=Q-D[m-1,:]-D[N-m,:]
      S1[m]=Q/(N-m)
  return np.mean(S1-2*S2,axis=1)

traj = gsdhoomd.HOOMDTrajectory(fl.GSDFile(sys.argv[1],'rb'))
box = traj[-1].configuration.box
fbox = Box(Lx=box[0],Ly=box[1],Lz=box[2],xy=box[3],xz=box[4],yz=box[5])
positions_unwrap = []
for f in traj:
    pos_unwrap = np.copy(f.particles.position)
    fbox.unwrap(pos_unwrap, f.particles.image)
    positions_unwrap.append(pos_unwrap)

msd = msd_fft_n(np.array(positions_unwrap))
with open('msd.dat','w') as f:
    for tau,val in enumerate(msd):
        f.write('{} {}\n'.format(tau,val))
