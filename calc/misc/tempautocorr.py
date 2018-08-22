## COMPUTE AND PRODUCE TIMESCALE PLOTS
from __future__ import print_function
path = '/home/mkloewer/python/swm/'
import os; os.chdir(path) # change working directory
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

def acf(x,l):
    """ autocorrelation function of vector x up to lag l."""
    return np.array([1]+[np.corrcoef(x[:-i],x[i:])[0,1] for i in range(1,l)])

def findc(x,a):
    """ find crossing of vector x with value a."""
    return np.argmin(abs(x-a))

runfolders = [7,8,9,10]

p1 = np.array([300,1920])*1e3
p2 = np.array([2880,1920])*1e3
pi = np.zeros((3,2,2),dtype=np.int)      # (u,v,T) x (p1,p2) x (i,j)

## read data
for r,i in zip(runfolders,range(len(runfolders))):
    runpath = path+'data/run%04i' % r
    param = np.load(runpath+'/param.npy').all()
    
    # find p1 and p2 indices
    for ig,g in enumerate(['u','v','T']):
        for ip,p in enumerate([p1,p2]):
            for ij,xy in enumerate(['x','y']):
                pi[ig,ip,ij] = findc(param[xy+'_'+g],p[ij])
        
    ncu = Dataset(runpath+'/u.nc')
    ncv = Dataset(runpath+'/v.nc')
    nch = Dataset(runpath+'/h.nc')
    
    istart = 0
    
    if i == 0:
        u = ncu.variables['u'][istart:,pi[0,:,0],pi[0,:,1]][:,[0,1],[0,1]]
        v = ncv.variables['v'][istart:,pi[1,:,0],pi[1,:,1]][:,[0,1],[0,1]]
        h = nch.variables['h'][istart:,pi[2,:,0],pi[2,:,1]][:,[0,1],[0,1]]
        t = nch.variables['t'][istart:]

    else:
        u = np.concatenate((u,ncu.variables['u'][1:,pi[0,:,0],pi[0,:,1]][:,[0,1],[0,1]]))
        v = np.concatenate((v,ncv.variables['v'][1:,pi[1,:,0],pi[1,:,1]][:,[0,1],[0,1]]))
        h = np.concatenate((h,nch.variables['h'][1:,pi[2,:,0],pi[2,:,1]][:,[0,1],[0,1]]))
        t = np.hstack((t,nch.variables['t'][1:]))
    
    ncu.close()
    ncv.close()
    nch.close()
    
    print('run %i read.' % r)

## computation
l = 200     # in 1/4 days

acfs = np.zeros((l,3,2))
for iv,var in enumerate([u,v,h]):
    for ip in range(2):
        acfs[:,iv,ip] = acf(var[:,ip],l)
        
dt = t[1]-t[0]
time = np.arange(l)*dt/24/3600

## STORING
dic = dict()
all_var2export = ['time','acfs','p1','p2']

for var in all_var2export:
    exec('dic[var] ='+var)
    
np.save(runpath+'/analysis/acfs.npy',dic)
print('Everything stored.')
