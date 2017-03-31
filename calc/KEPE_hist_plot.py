## HISTOGRAM PLOTTING FOR KINETIC AND POTENTIAL ENERGY
from __future__ import print_function
path = '/home/mkloewer/python/swm/'
import os; os.chdir(path) # change working directory
import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt
import time as tictoc
from netCDF4 import Dataset
import glob
from matplotlib.colors import BoundaryNorm,LogNorm

# OPTIONS
runfolder = 3
print(('Creating KE PE histogram plot from run %i') % runfolder)

## read data
runpath = path+'data/run%04i' % runfolder

D = np.load(runpath+'/analysis/KEPE_hist.npy').all()
for k in list(D.keys()):
   exec(k+' = D["'+k+'"]')
   
print('NPY data read.')

## PLOT

fig,ax1 = plt.subplots(1,1)

ax1.contour(D['KE_mid'],D['PE_mid'],D['KEPEH'],20,cmap='viridis')
plt.show()
