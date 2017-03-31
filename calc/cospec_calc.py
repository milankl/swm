## COMPUTE COSPEC
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
print('Calculating cospecs from run ' + str(runfolder))
    
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
        #h = np.load(runpath+'/h_sub.npy')        
        time = np.load(runpath+'/t_sub.npy')
        print('run %i read.' % r)

    else:
        u = np.concatenate((u,np.load(runpath+'/u_sub.npy')))
        v = np.concatenate((v,np.load(runpath+'/v_sub.npy')))
        #h = np.concatenate((h,np.load(runpath+'/h_sub.npy')))
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

lats = range(10,125,5)
lats = [64]
#lons = [64]
for lat in lats:
#for lon in lons:
    #f,k,p = cospec(u[:,:,lon],param['dy'],dt)
    f,k,pu = cospec(u[:,lat,:],param['dx'],dt)
    f,k,pv = cospec(v[:,lat,:-1],param['dx'],dt)
    
    p_all.append(pu+pv)
    print(lat)
    #print(lon)

p = np.array(p_all).mean(axis=0)

##

dic = dict()
all_var2export = ['f','k','p']

for v in all_var2export:
    exec('dic[v] ='+v)
    
np.save(runpath+'/analysis/cospec_eke_1lat.npy',dic)
print('Everything stored.')