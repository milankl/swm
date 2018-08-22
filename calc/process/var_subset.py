## READ VARIABLE FROM SEVERAL NCFILES and store subset of it as NPY
from __future__ import print_function

path = '/home/mkloewer/python/swm/'

import os; os.chdir(path) # change working directory
import numpy as np
from netCDF4 import Dataset


# OPTIONS
runfolder = [14]
sub = 4

for r in runfolder:
    print(('Subsampling from run %i') % r)

    # read data
    runpath = path+'data/run%04i' % r

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
    time = nceta['t'][::sub]   # in seconds
    t = time / 3600. / 24.  # in days
    nceta.close()
    print('eta read.')
    np.save(runpath+'/eta_sub.npy',eta)
    np.save(runpath+'/t_sub.npy',time)
    del eta
    del time,t
