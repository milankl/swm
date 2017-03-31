## PRODUCE ENERGY CYCLE PLOTS
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

KEl = []
EKEl = []
MKEl = []

PEl = []
EPEl = []
MPEl = []

atime = []

R = 20  # reynolds averaging period in days

## functions

def rmean(x,l):
    return np.convolve(x,np.ones(l)/l)


## OPTIONS
for runfolder in [5,6]:

    ## read data
    runpath = path+'data/run%04i' % runfolder
    ncu = Dataset(runpath+'/u.nc')
    ncv = Dataset(runpath+'/v.nc')
    nch = Dataset(runpath+'/h.nc')
    
    u = ncu['u'][:][:,:,:]
    v = ncv['v'][:][:,:,:]
    h = nch['h'][:][:,:,:]
    time = ncu['t'][:][:]   # in seconds
    t = time / 3600. / 24.  # in days
    print('netCDF data read.')
    
    # close netcdfs
    ncu.close()
    ncv.close()
    nch.close()
    
    # read param
    global param
    param = np.load(runpath+'/param.npy').all()
    param['dat_type'] = np.float32
    
    # import functions
    exec(open(path+'swm_param.py').read())
    exec(open(path+'swm_operators.py').read())
    exec(open(path+'swm_output.py').read())
    param['output'] = 0
    
    set_grad_mat()
    set_interp_mat()
    set_lapl_mat()
    set_coriolis()
    set_forcing()
    
    tlen = len(time)
    Rp = int(R*3600.*24 / param['output_dt'])   # reynolds period
    Rp += (Rp+1) % 2    # make sure Rp is odd
    filt = np.ones(Rp) / float(Rp)
    mask = np.hstack((np.ones(Rp//2),np.zeros(tlen-Rp+1),np.ones(Rp//2)))
    
    ## create ouputfolder
    try:
        os.mkdir(runpath+'/plots')
    except:
        pass
    
    ## reshape u,v
    u = u.reshape((tlen,param['Nu'])).T
    v = v.reshape((tlen,param['Nv'])).T
    h = h.reshape((tlen,param['NT'])).T
    print('Reshape done.')
    
    ##
    PEm = (.5*param['g']*param['rho']*h**2).mean(axis=0)
    hbar = np.empty_like(h)
    for i in range(param['NT']):
        hbar[i,:] = np.convolve(h[i,:],filt,mode='same')
    
    MPE = (.5*param['g']*param['rho']*np.ma.masked_array(hbar,mask=np.array([mask]*param['NT']))**2).mean(axis=0)
    EPE = np.convolve(PEm,filt,mode='same') - MPE    
    
    KEm = (.5*param['rho']*param['H']*(IuT.dot(u**2) + IvT.dot(v**2))).mean(axis=0)
    del v
    print('Kinetic Energy done.')
    
    PEl.append(PEm)
    MPEl.append(MPE)
    EPEl.append(EPE)
    KEl.append(KEm)
    atime.append(t)

##

fig,ax1 = plt.subplots(1,1,sharex=True,figsize=(8,8))

colors = ['b','b','g','g','r','r']
labels = [None,'dx=7.5km',None,'dx=30km',None,'dx=60km']

# do not show initial conditions
for i in range(len(KEl)):
    ax1.plot(KEl[i],PEl[i],colors[i],label=labels[i],lw=4)
    ax1.plot(KEl[i],MPEl[i],colors[i],label=labels[i],lw=2)
    ax1.plot(KEl[i],EPEl[i],colors[i])
    
ax1.legend(loc='best')
ax1.set_xlabel(r'KE [Jm$^{-2}$]')
ax1.set_ylabel(r'PE [Jm$^{-2}$]')

plt.tight_layout()
plt.show()