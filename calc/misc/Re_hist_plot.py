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

plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'

# OPTIONS
# runfolder = [0,6,10,14,15]
# print('Creating ReRo histogram plot from run '+str(runfolder)) 
# 
# # read data
# runpath = path+'data/run%04i' % runfolder[0]
# D1 = np.load(runpath+'/analysis/Re_hist.npy').all()   
# 
# runpath = path+'data/run%04i' % runfolder[1]
# D2 = np.load(runpath+'/analysis/Re_hist.npy').all()
# 
# runpath = path+'data/run%04i' % runfolder[2]
# D3 = np.load(runpath+'/analysis/Re_hist.npy').all()
# 
# runpath = path+'data/run%04i' % runfolder[3]
# D4 = np.load(runpath+'/analysis/Re_hist.npy').all()
# 
# runpath = path+'data/run%04i' % runfolder[4]
# D5 = np.load(runpath+'/analysis/Re_hist.npy').all()


# OPTIONS
runfolder = [3,10,13,12,14]
print('Creating ReRo histogram plot from run '+str(runfolder)) 

## read data
runpath = path+'data/newold/run%04i' % runfolder[0]
D1 = np.load(runpath+'/analysis/Re_hist.npy').all()   

runpath = path+'data/newold/run%04i' % runfolder[1]
D2 = np.load(runpath+'/analysis/Re_hist.npy').all()

runpath = path+'stoch/data/run%04i' % runfolder[2]
D3 = np.load(runpath+'/analysis/Re_hist.npy').all()

runpath = path+'stoch/data/run%04i' % runfolder[3]
D4 = np.load(runpath+'/analysis/Re_hist.npy').all()

runpath = path+'stoch/data/run%04i' % runfolder[4]
D5 = np.load(runpath+'/analysis/Re_hist.npy').all()

## PLOT
fig,ax1 = plt.subplots(1,1,figsize=(8,6))

ax1.plot(D1['Re_mid'],D1['ReH'],'C0',label=r'Low resolution, $\Delta x = $30km',lw=3)
ax1.plot(D2['Re_mid'],D2['ReH']/16,'C2',label=r'High resolution, $\Delta x = $7.5km',lw=3)
ax1.plot(D3['Re_mid'],D3['ReH'],'C3',label=r'LR + weak backscatter',ls='--')
ax1.plot(D4['Re_mid'],D4['ReH'],'C1',label=r'LR + moderate backscatter',ls='--')
ax1.plot(D5['Re_mid'],D5['ReH'],'C5',label=r'LR + strong backscatter',ls='--')

ax1.axvline(np.log10(D1['Re_mean']),c='C0',ls='-',lw=2)
ax1.axvline(np.log10(D2['Re_mean']),c='C2',ls='-',lw=2)
ax1.axvline(np.log10(D3['Re_mean']),c='C3',ls='--')
ax1.axvline(np.log10(D4['Re_mean']),c='C1',ls='--')
ax1.axvline(np.log10(D5['Re_mean']),c='C5',ls='--')

ax1.text(np.log10(D1['Re_mean']),5e5,'mean($R_e$)',rotation=90,ha='right',color='k')

#ax1.set_yscale('log')

ax1.legend(loc=2)

ax1.set_xlim(-3,5)
ax1.set_ylim(1,3e6)

ax1.set_title('Reynolds number histogram',loc='left')
ax1.set_xlabel('log$_{10}(R_e)$')
ax1.set_ylabel(r'$N$')

plt.tight_layout()
plt.savefig(path+'compare/Re_hist_nobf.png')
plt.close(fig)