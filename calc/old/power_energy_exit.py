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

Exu = []    # power dissipation u
Exv = []    # and v

Exui = []   # energy dissipation u
Exvi = []   # and v

atime = []

## OPTIONS
for runfolder in [1,2,3,4,5,6]:

    ## read data
    runpath = path+'data/run%04i' % runfolder
    ncu = Dataset(runpath+'/u.nc')
    ncv = Dataset(runpath+'/v.nc')  
    u = ncu['u'][:][:,:,:]
    v = ncv['v'][:][:,:,:]
    time = ncu['t'][:][:]   # in seconds
    t = time / 3600. / 24.  # in days
    print('netCDF data read.')
    
    # close netcdfs
    ncu.close()
    ncv.close()
    
    # read param
    global param
    param = np.load(runpath+'/param.npy').all()
    param['dat_type'] = np.float64
    
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
    ## create ouputfolder
    try:
        os.mkdir(runpath+'/plots')
    except:
        pass
    
    ## reshape u,v
    u = u.reshape((tlen,param['Nu'])).T
    v = v.reshape((tlen,param['Nv'])).T
    print('Reshape done.')

    diss_u_m = (param['B']*(u*LLu.dot(u))*param['rho']*param['H']).mean(axis=0)
    diss_uim = np.cumsum(diss_u_m*param['output_dt'],axis=0)
    
    diss_v_m = (param['B']*(v*LLv.dot(v))*param['rho']*param['H']).mean(axis=0)
    diss_vim = np.cumsum(diss_v_m*param['output_dt'],axis=0)
    
    del u,v
    print('Exit Energy done.')
    
    Exu.append(diss_u_m)
    Exui.append(diss_uim)
    
    Exv.append(diss_v_m)
    Exvi.append(diss_vim)
    atime.append(t)

##

Exui[1] = Exui[1] + Exui[0][-1]
Exui[3] = Exui[3] + Exui[2][-1]
Exui[5] = Exui[5] + Exui[4][-1]

Exvi[1] = Exvi[1] + Exvi[0][-1]
Exvi[3] = Exvi[3] + Exvi[2][-1]
Exvi[5] = Exvi[5] + Exvi[4][-1]

##

fig,(ax1,ax2) = plt.subplots(2,1,sharex=True,figsize=(8,8))

colors = ['b','b','g','g','r','r']
labels = [None,'dx=7.5km',None,'dx=30km',None,'dx=60km']

for i in range(len(Exu)):
    ax1.plot(atime[i],Exu[i],colors[i],label=labels[i])
    ax1.plot(atime[i],Exv[i],colors[i]+'--')
    
    ax2.plot(atime[i],Exui[i],colors[i],label=labels[i])
    ax2.plot(atime[i],Exvi[i],colors[i]+'--')    

ax1.plot([atime[0][0],atime[1][-1]],[0,0],'grey')
ax1.legend(loc='best')
ax2.set_xlabel(r'time [days]')
ax1.set_ylabel(r'Power dissipation [Wm$^{-2}$]')
ax2.set_ylabel(r'Energy dissipation [Jm$^{-2}$]')

plt.tight_layout()
plt.show()