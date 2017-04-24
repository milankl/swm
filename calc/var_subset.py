## READ VARIABLE FROM SEVERAL NCFILES and store subset of it as NPY
from __future__ import print_function

# path
import os
path = os.path.dirname(os.getcwd()) + '/'   # on level above
os.chdir(path)                              # change working directory

import numpy as np
from netCDF4 import Dataset

# OPTIONS
runfolder = [1,2]

for r in runfolder:    
    print(('Subsampling from run %i') % r)
    
    runpath = path+'data/run%04i' % r
    sub = 4     # read only every sub-th time step 
    
    ncu = Dataset(runpath+'/u.nc')
    u = ncu['u'][:][::sub,:,:]
    ncu.close()
    print('u read.')
    np.save(runpath+'/u_sub.npy',u)
    del u
    
    ncv = Dataset(runpath+'/v.nc')
    v = ncv['v'][:][::sub,:,:]
    ncv.close()
    print('v read.')
    np.save(runpath+'/v_sub.npy',v)
    del v
    
    nceta = Dataset(runpath+'/eta.nc')
    eta = nceta['eta'][:][::sub,:,:]
    t = nceta['t'][::sub]   # in seconds
    nceta.close()
    print('eta read.')
    np.save(runpath+'/eta_sub.npy',h)
    np.save(runpath+'/t_sub.npy',t)
    del h,t