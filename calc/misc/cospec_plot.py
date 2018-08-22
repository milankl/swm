## COSPEC PLOT
from __future__ import print_function
path = '/home/mkloewer/python/swm/'
import os; os.chdir(path) # change working directory
import numpy as np
from scipy import sparse
import time as tictoc
from netCDF4 import Dataset
import glob
import matplotlib.pyplot as plt
from cmocean import cm

# OPTIONS
runfolder = [3,10]
var = 'eke'
print('Compare cospec plots from run ' + str(runfolder))

#lat1 = '_1lat'
lat1 = ''

runpath1 = path+'data/run%04i' % runfolder[0]
D1 = np.load(runpath1+'/analysis/cospec_'+var+lat1+'.npy').all()
param1 = np.load(runpath1+'/param.npy').all()

runpath2 = path+'data/run%04i' % runfolder[1]
D2 = np.load(runpath2+'/analysis/cospec_'+var+lat1+'.npy').all()
param2 = np.load(runpath2+'/param.npy').all()

##
Ro_max = param1['c_phase']/(param1['f_0'] - param1['beta']*param1['Ly']/2.)
Ro_min = param1['c_phase']/(param1['f_0'] + param1['beta']*param1['Ly']/2.)

# dispersion relation, divide omega and k by 2pi
wrossby_max = lambda k: -param1['beta']*k/(2*np.pi)/(2*(k/(2*np.pi))**2 + 1./(Ro_max)**2)/2./np.pi
wrossby_min = lambda k: -param1['beta']*k/(2*np.pi)/(2*(k/(2*np.pi))**2 + 1./(Ro_min)**2)/2./np.pi

## PLOTTING
sd = 3600.*24.  # from seconds to days for frequency axis
kf = 1e3        # wavenumber factor for plotting   

fig,(ax1,ax2) = plt.subplots(1,2,sharex=True,sharey=True,figsize=(14,7))
c1 = ax1.contourf(D1['k']*kf,D1['f']*sd,np.log10(D1['p']),np.linspace(4,10,51),cmap=cm.thermal,extend='both')
c2 = ax2.contourf(D2['k']*kf,D2['f']*sd,np.log10(D2['p']),np.linspace(4,10,51),cmap=cm.thermal,extend='both')

for D,ax,greyvalue in zip([D1,D2],[ax1,ax2],['0.7','0.3']):
    ax.plot(np.zeros(2),[D['f'][0]*sd,D['f'][-1]*sd],'w',lw=22)

    ax.plot(D['k']*kf,wrossby_min(D['k'])*sd,'c',label=r'$\frac{1}{2\pi}\omega_{Ro,min}(\frac{k}{2\pi})$',lw=1.5)
    ax.plot(D['k']*kf,wrossby_max(D['k'])*sd,'c--',label=r'$\frac{1}{2\pi}\omega_{Ro,max}(\frac{k}{2\pi})$',lw=1.5)
    
    ax.plot(1./Ro_min*kf*np.ones(2),[D['f'][0]*sd,D['f'][-1]*sd],greyvalue,label=r'$Ro_{min}$',lw=1.5)
    ax.plot(1./Ro_max*kf*np.ones(2),[D2['f'][0]*sd,D['f'][-1]*sd],greyvalue,label=r'$Ro_{max}$',linestyle='dashed',lw=1.5)
    
    ax.plot(-1./Ro_min*kf*np.ones(2),[D['f'][0]*sd,D['f'][-1]*sd],greyvalue,lw=1.5)
    ax.plot(-1./Ro_max*kf*np.ones(2),[D2['f'][0]*sd,D['f'][-1]*sd],greyvalue,linestyle='dashed',lw=1.5)
    
    ytick = np.array([200,50,20,15,10,8,7,6,5,4,3])
    ax.set_yticks(1./ytick)
    ax.set_yticklabels(ytick)
    
    xtick = np.array([-2,-5,-10,-40,40,10,5,2])*100
    ax.set_xticks(1./xtick)
    ax.set_xticklabels(xtick)
    
    ax.tick_params(axis='x', color='white')
    ax.tick_params(axis='y', color='white')
    
    ax.set_xticklabels(ax.xaxis.get_majorticklabels(), rotation=90)

ax1.set_xlim(-0.005,0.005)
ax1.set_ylim(0,.34)

ax1.legend(loc=1)

ax1.set_xlabel(r'wavelength [km]')
ax2.set_xlabel(r'wavelength [km]')
ax1.set_ylabel('period [days]')

if var == 'h':
    varname = '$\eta$'
elif var == 'u':
    varname = '$u$'
elif var == 'v':
    varname = '$v$'
elif var == 'eke':
    varname = 'EKE'

ax1.set_title('2D power spectrum of '+varname+': Low resolution $\Delta x = 30$km')
ax2.set_title('High resolution $\Delta x = 7.5$km')
plt.tight_layout()
plt.savefig(path+'/compare/cospec_'+var+lat1+'.png')
plt.close(fig)


