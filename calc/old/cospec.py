## COMPUTE TRISPEC
from __future__ import print_function
path = '/home/mkloewer/python/swm/'
import os; os.chdir(path) # change working directory
import numpy as np
from scipy import sparse
import time as tictoc
from netCDF4 import Dataset
import glob
import matplotlib.pyplot as plt

# OPTIONS
runfolder = [2,3]
print('Calculating histograms from run ' + str(runfolder))
    
def cospec(A,dx,dt):
    """ assumed shape of A is time x space. """
    nt,nx = np.shape(A)
    
    f = np.fft.fftfreq(nt,dt)
    k = np.fft.fftfreq(nx,dx)
    
    idxf = int(np.ceil(nt/2.))    # index of maximum positive frequency
    idxk = int(np.ceil(nx/2.))    # index of maximum positive wavenumber

    # 2D FFT
    p = abs(np.fft.fft2(A))**2
    
    f = f[1:idxf]    # kill 0 and negative frequencies
    # reorder to have -kmax to +kmax
    k = np.hstack((k[idxk:],k[1:idxk]))
    # apply also for p
    p = np.hstack((p[1:idxf,idxk:],p[1:idxf,1:idxk]))[:,::-1]
    
    #TODO issue: p_max is shifted by -1 in positive k-direction
    return f,k,p

## read data
for r,i in zip(runfolder,range(len(runfolder))):
    runpath = path+'data/run%04i' % r
    
    if i == 0:
        u = np.load(runpath+'/u_sub.npy')
        v = np.load(runpath+'/v_sub.npy')
        h = np.load(runpath+'/h_sub.npy')        
        time = np.load(runpath+'/t_sub.npy')
        print('run %i read.' % r)

    else:
        u = np.concatenate((u,np.load(runpath+'/u_sub.npy')))
        v = np.concatenate((v,np.load(runpath+'/v_sub.npy')))
        h = np.concatenate((h,np.load(runpath+'/h_sub.npy')))
        time = np.hstack((time,np.load(runpath+'/t_sub.npy')))
        print('run %i read.' % r)

t = time / 3600. / 24.  # in days
dt = time[1] - time[0]

## read param
global param
param = np.load(runpath+'/param.npy').all()

# import functions
exec(open(path+'swm_param.py').read())
exec(open(path+'swm_operators.py').read())
exec(open(path+'swm_output.py').read())
param['output'] = 0

set_grad_mat()
set_interp_mat()
set_lapl_mat()
set_coriolis()

tlen = len(time)

##
p_all = []

lats = range(10,110,5)
for lat in lats:
    f,k,p = cospec(h[:,lat,:],param['dx'],dt)
    p_all.append(p)
    print(lat)

p = np.array(p_all).mean(axis=0)

# rossby radius
Ro_max = param['c_phase']/h2mat(f_T)[lats[0],0]
Ro_min = param['c_phase']/h2mat(f_T)[lats[-1],0]

# dispersion relation
wrossby_max = -param['beta']*k/(k**2 + 1./Ro_max**2)
wrossby_min = -param['beta']*k/(k**2 + 1./Ro_min**2)

##
# plotting 2D
sd = 3600.*24.  # from seconds to days for frequency axis
kf = 1e5        # wavenumber factor for plotting   

fig,ax = plt.subplots(1,1)
c1 = ax.contourf(k*kf,f*sd,np.log10(p),np.linspace(0,10,51),cmap='viridis',extend='both')
plt.colorbar(c1,ax=ax)

ax.plot(k*kf,wrossby_min*sd,'k')
ax.plot(k*kf,wrossby_max*sd,'k')

ax.plot(1./Ro_min*kf*np.ones(2),[f[0]*sd,f[-1]*sd],'w')
ax.plot(1./Ro_max*kf*np.ones(2),[f[0]*sd,f[-1]*sd],'w--')

ax.plot(-1./Ro_min*kf*np.ones(2),[f[0]*sd,f[-1]*sd],'w')
ax.plot(-1./Ro_max*kf*np.ones(2),[f[0]*sd,f[-1]*sd],'w--')

ax.set_xlim(-0.5,0.5)
ax.set_ylim(0,.2)

ytick = np.array([200,50,20,10,5])
ax.set_yticks(1./ytick)
ax.set_yticklabels(ytick)

xtick = np.array([-2,-5,-10,-40,40,10,5,2])
ax.set_xticks(1./xtick)
ax.set_xticklabels(xtick)

ax.set_xlabel(r'wavelength [100km]')
ax.set_ylabel('period [days]')
ax.set_title('2D power spectrum of $\eta$')
plt.tight_layout()
plt.savefig(runpath+'/plots/cospec/cospec_h.png')
plt.close(fig)


