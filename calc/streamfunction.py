# PLOT STREAMFUNCTION
from __future__ import print_function
path = '/home/mkloewer/python/swm/'
import os; os.chdir(path) # change working directory
import numpy as np
from scipy import sparse
from scipy.integrate import cumtrapz
import matplotlib.pyplot as plt
import time as tictoc
from netCDF4 import Dataset
import glob
from cmocean import cm

## OPTIONS
runfolder = [3,10]
print('Compare mean plots from run ' + str(runfolder))

## read data

runpath1 = path+'data/run%04i' % runfolder[0]
D1 = np.load(runpath1+'/analysis/mean.npy').all()
param1 = np.load(runpath1+'/param.npy').all()

# runpath2 = path+'data/run%04i' % runfolder[1]
# D2 = np.load(runpath2+'/analysis/mean.npy').all()
# param2 = np.load(runpath2+'/param.npy').all()
# 
# runpath3 = path+'stoch/data/run%04i' % 5
# D3 = np.load(runpath3+'/analysis/mean.npy').all()
# param3 = np.load(runpath3+'/param.npy').all()
    
## CALCULATE STREAMFUNCTION
global param
param = np.load(runpath1+'/param.npy').all()
exec(open(path+'swm_operators.py').read())
exec(open(path+'swm_output.py').read())
param['output'] = 0
set_grad_mat()
set_interp_mat()

PSI1v = cumtrapz(v2mat(D1['vm']),axis=1,dx=param['dx'],initial=0)
PSI1u = -cumtrapz(u2mat(D1['um']),axis=0,dx=param['dy'],initial=0)

levs = np.linspace(-abs(PSI1u).max(),abs(PSI1u).max(),64)

fig,(ax1,ax2) = plt.subplots(1,2,sharex=True,sharey=True)

ax1.contourf(param['x_v'],param['y_v'],PSI1v,levs,cmap=cm.balance)
ax2.contourf(param['x_u'],param['y_u'],PSI1u,levs,cmap=cm.balance)

plt.show()