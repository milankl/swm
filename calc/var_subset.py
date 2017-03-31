## READ VARIABLE FROM SEVERAL NCFILES and store subset of it as NPY
from __future__ import print_function
path = '/home/mkloewer/github/swm/'
import os; os.chdir(path) # change working directory
import numpy as np
from scipy import sparse
import time as tictoc
from netCDF4 import Dataset
import glob

# OPTIONS
runfolder = [0,1,2]

for r in runfolder:    
    print(('Subsampling from run %i') % r)
    
    ## read data
    runpath = path+'data/run%04i' % r
    
    ##
    sub = 4
    
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
    
    nch = Dataset(runpath+'/h.nc')
    h = nch['h'][:][::sub,:,:]
    time = nch['t'][::sub]   # in seconds
    t = time / 3600. / 24.  # in days
    nch.close()
    print('h read.')
    np.save(runpath+'/h_sub.npy',h)
    np.save(runpath+'/t_sub.npy',time)
    del h
    del time,t