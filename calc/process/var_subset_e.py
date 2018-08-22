## READ VARIABLE FROM SEVERAL NCFILES and store subset of it as NPY
from __future__ import print_function

path = '/home/mkloewer/python/swm/'

import os; os.chdir(path) # change working directory
import numpy as np
from netCDF4 import Dataset

# OPTIONS
runfolder = [10,11,12,13,14,15]
sub = 4

for r in runfolder:
    print(('Subsampling from run %i') % r)

    ## read data
    runpath = path+'data/run%04i' % r
    nce = Dataset(runpath+'/e.nc')
    e = nce['e'][:][::sub,:,:]
    nce.close()

    # save subset as .npy
    np.save(runpath+'/e_sub.npy',e)
    del e
