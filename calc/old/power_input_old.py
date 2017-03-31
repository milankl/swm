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

FElist = []
FEilist = []
atime = []

## OPTIONS
for runfolder in [1,2,3,4,5,6]:

    ## read data
    runpath = path+'data/run%04i' % runfolder
    ncu = Dataset(runpath+'/u.nc')    
    u = ncu['u'][:][:,:,:]
    time = ncu['t'][:][:]   # in seconds
    t = time / 3600. / 24.  # in days
    print('netCDF data read.')
    
    # close netcdfs
    ncu.close()
    
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
    print('Reshape done.')

    FEm = ((u.T*Fx).T*param['rho']*param['H']).mean(axis=0)
    FEim = np.cumsum((u.T*Fx).T*param['rho']*param['H']*param['output_dt'],axis=1).mean(axis=0)
    del u
    print('Input Energy done.')
    
    FElist.append(FEm)
    FEilist.append(FEim)
    atime.append(t)

##

FEilist[1] = FEilist[1] + FEilist[0][-1]
FEilist[3] = FEilist[3] + FEilist[2][-1]
FEilist[5] = FEilist[5] + FEilist[4][-1]

##

fig,(ax1,ax2) = plt.subplots(2,1,sharex=True,figsize=(8,8))

colors = ['b','b','g','g','r','r']
labels = [None,'dx=7.5km',None,'dx=30km',None,'dx=60km']

for i in range(len(FElist)):
    ax1.plot(atime[i],FElist[i],colors[i],label=labels[i])
    ax2.plot(atime[i],FEilist[i],colors[i],label=labels[i])
    

ax1.plot([atime[0][0],atime[1][-1]],[0,0],'grey')
ax1.legend(loc='best')
ax2.set_xlabel(r'time [days]')
ax1.set_ylabel(r'Power input [Wm$^{-2}$]')
ax2.set_ylabel(r'Energy input [Jm$^{-2}$]')

plt.tight_layout()
plt.show()