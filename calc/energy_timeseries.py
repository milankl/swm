## PLOT ENERGY AND ENSTROPY TIME SERIES
from __future__ import print_function
path = '/home/mkloewer/python/swm/'
import os; os.chdir(path) # change working directory
import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt
import time as tictoc
from netCDF4 import Dataset
import glob
import cmocean
import matplotlib.patches as patches

# OPTIONS
runfolder = [1,3,5,6,10]
print('Compare mean plots from run ' + str(runfolder))

## read data

runpath1 = path+'data/run%04i' % runfolder[0]
D1 = np.load(runpath1+'/analysis/mean.npy').all()

runpath2 = path+'data/run%04i' % runfolder[1]
D2 = np.load(runpath2+'/analysis/mean.npy').all()

KEm_lo = np.hstack((D1['KEm'],D2['KEm'])) * 1e-3
PEm_lo = np.hstack((D1['PEm'],D2['PEm'])) * 1e-3
Zm_lo = np.hstack((D1['Zm'],D2['Zm'])) * 1e11
time_lo = np.hstack((D1['time'],D2['time'])) / 3600. / 24. / 365. # in years

del D1,D2

runpath3 = path+'data/run%04i' % runfolder[2]
D3 = np.load(runpath3+'/analysis/mean.npy').all()

runpath4 = path+'data/run%04i' % runfolder[3]
D4 = np.load(runpath4+'/analysis/mean.npy').all()

runpath5 = path+'data/run%04i' % runfolder[4]
D5 = np.load(runpath5+'/analysis/mean.npy').all()

KEm_hi = np.hstack((D3['KEm'],D4['KEm'],D5['KEm'])) * 1e-3
PEm_hi = np.hstack((D3['PEm'],D4['PEm'],D5['PEm'])) * 1e-3
Zm_hi = np.hstack((D3['Zm'],D4['Zm'],D5['Zm'])) * 1e11
time_hi = np.hstack((D3['time'],D4['time'],D5['time'])) / 3600. / 24. / 365. # in years

del D3,D4,D5

## PLOTTING   

fig,(ax1,ax2,ax3) = plt.subplots(3,1,sharex=True)

ax1.plot(time_hi,KEm_hi,'g',label=r'High resolution, $\Delta x = 7.5$km')
ax1.plot(time_lo,KEm_lo,'b',label=r'Low resolution, $\Delta x = 30$km')


ax2.plot(time_lo,PEm_lo)
ax2.plot(time_hi,PEm_hi)

ax3.plot(time_lo,Zm_lo)
ax3.plot(time_hi,Zm_hi)

ax1.add_patch(patches.Rectangle((0,0),10,ax1.get_ylim()[1],color='grey',alpha=.3))
ax2.add_patch(patches.Rectangle((0,0),10,ax1.get_ylim()[1],color='grey',alpha=.3))
ax3.add_patch(patches.Rectangle((0,0),10,ax1.get_ylim()[1],color='grey',alpha=.3))
#invisible bar for legend
ax1.bar(0,0,color='grey',alpha=.3,label='spin-up')

ax1.set_title('Kinetic Energy')
ax1.set_ylabel(r'kJm$^{-2}$')

ax2.set_title('Potential Energy')
ax2.set_ylabel(r'kJm$^{-2}$')

ax3.set_title('Potential Enstrophy')
ax3.set_ylabel(r'$10^{-11}$m$^{-1}$s$^{-2}$')
ax1.legend(loc=5,fontsize=10,ncol=2)

plt.tight_layout()
#plt.savefig('compare/Energy_timeseries.png')
plt.show()