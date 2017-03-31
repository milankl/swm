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

from scipy.integrate import cumtrapz

## read data
runpath1 = path+'data/run%04i' % 3
runpath11 = path+'data/run%04i' % 1
D1 = np.load(runpath1+'/analysis/power_energy_timeseries.npy').all()
D1m = np.load(runpath11+'/analysis/mean.npy').all()
param1 = np.load(runpath1+'/param.npy').all()

runpath2 = path+'data/run%04i' % 10
D2 = np.load(runpath2+'/analysis/power_energy_timeseries.npy').all()

runpath3 = path+'data/run%04i' % 5
D3 = np.load(runpath3+'/analysis/mean.npy').all()

runpath4 = path+'data/run%04i' % 6
D4 = np.load(runpath4+'/analysis/mean.npy').all()

KEm_hi = np.hstack((D3['KEm'],D4['KEm']))
PEm_hi = np.hstack((D3['PEm'],D4['PEm']))
time_hi = np.hstack((D3['time'],D4['time']))

del D3,D4

## Plotting
s = 1./(24*3600*365)

fig,(ax1,ax2) = plt.subplots(2,1,sharex=True)

ax1.plot(D1['time']*s,D1['InPower'],'C0')
ax1.plot(D2['time']*s,D2['InPower'],'C2')

ax2.plot(D1['time']*s,D1['ExPower_u']+D1['ExPower_v'],'C0',label=r'Low resolution, $\Delta x = 30$km')
ax2.plot(D2['time']*s,D2['ExPower_u']+D2['ExPower_v'],'C2',label=r'High resolution, $\Delta x = 7.5$km')

ax1.set_xlim(0,5)
ax1.set_ylim(-0.02,0.04)
ax2.set_ylim(-0.02,0)
ax1.plot(ax1.get_xlim(),[0,0],'0.5')
ax2.plot(ax2.get_xlim(),[0,0],'0.5')

ax1.set_title('Power input by forcing')
ax1.set_ylabel(r'Power [Wm$^{-2}$]')

ax2.set_title('Dissipation power by mixing')
ax2.set_ylabel(r'Power [Wm$^{-2}$]')
ax2.set_xlabel('time [years]')

ax2.legend(loc=3,ncol=2)
plt.tight_layout()

n = D1m['time'].shape[0]
n2 = time_hi.shape[0]
en_res = (D1['InEnergy']+D1['ExEnergy'])[:n] - D1m['KEm'] - D1m['PEm']
en_res_hi = (D2['InEnergy']+D2['ExEnergy'])[:n2] - KEm_hi - PEm_hi

pow_res = (D1['InPower']+D1['ExPower_u']+D1['ExPower_v'])[:n] - np.gradient(D1m['KEm']+D1m['PEm'],np.gradient(D1m['time']))
pow_res_hi = (D2['InPower']+D2['ExPower_u']+D2['ExPower_v'])[:n2] - np.gradient(KEm_hi+PEm_hi,np.gradient(time_hi))

fig2,(ax21,ax22) = plt.subplots(2,1,sharex=True)


ax21.plot(D1m['time']*s,pow_res/np.gradient(D1m['KEm']+D1m['PEm'],np.gradient(D1m['time'])),'C0')
ax21.plot(time_hi*s,pow_res_hi/np.gradient(KEm_hi+PEm_hi,np.gradient(time_hi)),'C2')

ax22.plot(D1m['time']*s,en_res/(D1m['KEm']+D1m['PEm'])*100,'C0',label=r'Low resolution, $\Delta x = 30$km')
ax22.plot(time_hi*s,en_res_hi/(KEm_hi+PEm_hi)*100,'C2',label=r'High resolution, $\Delta x = 7.5$km')
ax21.set_xlim(0,10)
ax22.legend(loc=3)

ax21.set_title(r'Power budget residual: In + Out - $\frac{\partial}{\partial t}$KE - $\frac{\partial}{\partial t}$PE')
ax22.set_title(r'Energy budget residual: $\int$(In + Out)dt - KE - PE')
ax22.set_xlabel('time [years]')
ax22.set_ylabel('%')
ax21.set_ylabel('%')
ax21.set_ylim(-1,1)
ax22.set_ylim(-5,5)


plt.tight_layout()
plt.show()

