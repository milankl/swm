## READ VARIABLE FROM SEVERAL NCFILES and store subset of it as NPY
from __future__ import print_function

path = '/network/aopp/cirrus/pred/kloewer/swm_bf_cntrl/data/'
#path = '/network/aopp/cirrus/pred/kloewer/swm_back_ronew/'

import os; os.chdir(path) # change working directory
import numpy as np
from netCDF4 import Dataset

# OPTIONS
runfolder = [0,6]
s = 40              # read s-th last time step


for r in runfolder:
    print(('Store last time step from run %i') % r)

    ## read data
    runpath = path+'run%04i' % r

    ncu = Dataset(runpath+'/u.nc')
    u = ncu['u'][-s,:,:]
    ncu.close()
    print('u read.')
    np.save(runpath+'/u_last.npy',u)
    del u

    ncv = Dataset(runpath+'/v.nc')
    v = ncv['v'][-s,:,:]
    ncv.close()
    print('v read.')
    np.save(runpath+'/v_last.npy',v)
    del v

    nceta = Dataset(runpath+'/eta.nc')
    eta = nceta['eta'][-s,:,:]
    #time = nceta['t'][::sub]   # in seconds
    #t = time / 3600. / 24.  # in days
    nceta.close()
    print('eta read.')
    np.save(runpath+'/eta_last.npy',eta)
    del eta
