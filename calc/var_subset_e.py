## READ VARIABLE FROM SEVERAL NCFILES and store subset of it as NPY
from __future__ import print_function
path = '/home/mkloewer/python/swm/'
import os; os.chdir(path) # change working directory
import numpy as np
from scipy import sparse
import time as tictoc
from netCDF4 import Dataset
import glob

# OPTIONS
runfolder = [12,13,14]

for r in runfolder:    
    print(('Subsampling from run %i') % r)
    runpath = path+'stoch/data/run%04i' % r
    
    sub = 4    
    nce = Dataset(runpath+'/e.nc')
    e = nce['e'][:][::sub,:,:]
    nce.close()
    print('e read.')
    np.save(runpath+'/e_sub.npy',e)
    del e