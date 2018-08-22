## HISTOGRAM PLOTTING FOR REYNOLDS AND ROSSBY NUMBERS
from __future__ import print_function
path = '/home/mkloewer/python/swm/'
import os; os.chdir(path) # change working directory
import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt
import time as tictoc
from netCDF4 import Dataset
import glob
import matplotlib

# OPTIONS
runfolder = [3,10]
print('Creating ReRo histogram plot from run '+str(runfolder)) 

## read data
runpath = path+'data/run%04i' % runfolder[0]
D1 = np.load(runpath+'/analysis/ReRo_hist.npy').all()   
D11 = np.load(runpath+'/analysis/mean.npy').all()  

runpath = path+'data/run%04i' % runfolder[1]
D2 = np.load(runpath+'/analysis/ReRo_hist.npy').all()
D21 = np.load(runpath+'/analysis/mean.npy').all()  

runpath = path+'stoch/data/run%04i' % 12
D3 = np.load(runpath+'/analysis/ReRo_hist.npy').all()
D31 = np.load(runpath+'/analysis/mean.npy').all()  

## PLOT
fig,ax1 = plt.subplots(1,1,figsize=(8,6))

levs = np.array([12000,40000,106000])

ax1.contour(D2['Re_mid'],D2['Ro_mid'],D2['ReRoH'].T,levs*16.,colors='C2',linewidths=np.linspace(1,5,3),alpha=0.8)
ax1.contour(D1['Re_mid'],D1['Ro_mid'],D1['ReRoH'].T,levs,colors='C0',linewidths=np.linspace(1,5,3),alpha=0.8)
ax1.contour(D3['Re_mid'],D3['Ro_mid'],D3['ReRoH'].T,levs,colors='C1',linewidths=np.linspace(1,5,3),alpha=0.8)

maxij = [0,]*3
maxij[0] = np.unravel_index(np.argmax(D1['ReRoH']),(200,200))
maxij[1] = np.unravel_index(np.argmax(D2['ReRoH']),(200,200))
maxij[2] = np.unravel_index(np.argmax(D3['ReRoH']),(200,200))

Remaxs = [D['Re_mid'][ij[0]] for D,ij in zip([D1,D2,D3],maxij)]
Romaxs = [D['Ro_mid'][ij[1]] for D,ij in zip([D1,D2,D3],maxij)]

Remeans = [np.log10(D['Rem'].mean()) for D in [D11,D21,D31]]
Romeans = [np.log10(D['Rom'].mean()) for D in [D11,D21,D31]]

ax1.scatter(Remaxs,Romaxs,80,c=['C0','C2','C1'],edgecolors='k',zorder=2,label='mode')
ax1.scatter(Remeans,Romeans,80,c=['C0','C2','C1'],edgecolors='k',zorder=2,label='mean',marker='v')

# for legend
ax1.plot([0,0],[0,0],'C2',label=r'High resolution, $\Delta x = 7.5$km')
ax1.plot([0,0],[0,0],'C0',label=r'Low resolution, $\Delta x = 30$km')
ax1.plot([0,0],[0,0],'C1',label=r'Low resolution + backscatter')

ax1.legend(loc=4,scatterpoints=3)

ax1.set_xlim(0,4)
ax1.set_ylim(-3,0)
ax1.set_xlabel(r'$\log_{10}(Re)$')
ax1.set_ylabel(r'$\log_{10}(Ro)$')
ax1.set_title('2D Histogram of Reynolds and Rossby number')

plt.tight_layout()
plt.savefig(path+'compare/ReRo_2dhist.png')
plt.close(fig)