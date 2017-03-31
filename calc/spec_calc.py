## COMPUTE EKE SPECTRUM
from __future__ import print_function
path = '/home/mkloewer/github/swm/'
import os; os.chdir(path) # change working directory
import numpy as np
from scipy import sparse
import time as tictoc
from netCDF4 import Dataset
import glob
import matplotlib.pyplot as plt

exec(open(path+'swm_param.py').read())
exec(open(path+'swm_operators.py').read())
exec(open(path+'swm_rhs.py').read())
exec(open(path+'swm_integration.py').read())
exec(open(path+'swm_output.py').read())

# OPTIONS
runfolder = [2]
print('Calculating EKE spectrum from run ' + str(runfolder))

##
def eke_spec_avg(u,v,dx,dy):
    """ Computes a wavenumber-frequency plot for 3D (t,x,y) data via radial (k = sqrt(kx**2 + ky**2)) integration. TODO: correct normalisation, so that the integral in normal space corresponds to the integral in Fourier space.
    """
    
    nt,ny,nx = np.shape(u)
    kx = (1/(dx))*np.hstack((np.arange(0,(nx+1)/2.),np.arange(-nx/2.+1,0)))/float(nx)
    ky = (1/(dy))*np.hstack((np.arange(0,(ny+1)/2.),np.arange(-ny/2.+1,0)))/float(ny)

    kxx,kyy = np.meshgrid(kx,ky)
    # radial distance from kx,ky = 0[-700:]
    kk = np.sqrt(kxx**2 + kyy**2) 

    if nx >= ny: #kill negative wavenumbers
        k  = kx[:int(nx/2)+1]
    else:
        k  = ky[:int(ny/2)+1]

    dk = k[1] - k[0]

    # 2D FFT average
    p_eke = np.empty((nt,ny,nx))
    nxy2 = nx**2*ny**2

    
    for i in range(nt):
        pu = abs(np.fft.fft2(u[i,:,:]))**2/nxy2
        pv = abs(np.fft.fft2(v[i,:,:]))**2/nxy2
        p_eke[i,:,:] = pu+pv
        if ((i+1)/nt*100 % 5) < (i/nt*100 % 5):
            print(str(int((i+1)/nt*100.))+'%')
    
    p_eke_avg = .5*p_eke.mean(axis=0)

    # create radial coordinates, associated with k[i]
    # WRONG nearest point interpolation to get points within the -.5,.5 annulus
    rcoords = []
    for i in range(len(k)):
        #rcoords.append(np.where((kk>(k[i]-.5*dk))*(kk<=(k[i]+.5*dk))))
        rcoords.append(np.where(kk<k[i]))

    # mulitply by dk to have the corresponding integral
    eke_spec = np.zeros(len(k))
    for i in range(len(k)):
        eke_spec[i] = np.sum(p_eke_avg[rcoords[i][0],rcoords[i][1]])
    
    eke_spec = np.diff(eke_spec) / dk
    k = (k[:-1] + k[1:])/2.

    return k,eke_spec  # eliminate zero wavenumber

## read data
for r,i in zip(runfolder,range(len(runfolder))):
    runpath = path+'data/run%04i' % r
    
    if i == 0:
        u = np.load(runpath+'/u_sub.npy')[365:,:,:]
        v = np.load(runpath+'/v_sub.npy')[365:,:,:]
        #h = np.load(runpath+'/h_sub.npy')        
        time = np.load(runpath+'/t_sub.npy')[365:]
        print('run %i read.' % r)

    else:
        u = np.concatenate((u,np.load(runpath+'/u_sub.npy')))
        v = np.concatenate((v,np.load(runpath+'/v_sub.npy')))
        #h = np.concatenate((h,np.load(runpath+'/h_sub.npy')))
        time = np.hstack((time,np.load(runpath+'/t_sub.npy')))
        print('run %i read.' % r)

t = time / 3600. / 24.  # in days
tlen = len(time)
dt = time[1] - time[0]

## create ouputfolder
try:
    os.mkdir(runpath+'/analysis')
except:
   pass

## read param
global param
param = np.load(runpath+'/param.npy').all()
param['output'] = 0

set_grad_mat()
set_interp_mat()

## part 1 cut off boundary by 3 grid cells
k,p = eke_spec_avg(u[:,2:-3,2:-2],v[:,2:-2,2:-3],param['dx'],param['dy'])

dic = dict()
all_var2export = ['k','p']

for vars in all_var2export:
   exec('dic[vars] ='+vars)
    
np.save(runpath+'/analysis/spec_eke.npy',dic)

## include kinematic boundary condition on one side, pad with zero
u = np.pad(u,((0,0),(0,0),(1,1)),'constant')
v = np.pad(v,((0,0),(1,1),(0,0)),'constant')

u = np.pad(u,((0,0),(0,1),(0,0)),'edge')
v = np.pad(v,((0,0),(0,0),(0,1)),'edge')

u[:,-1,:] = -u[:,-1,:]
v[:,:,-1] = -v[:,:,-1]

k,p = eke_spec_avg(u,v,param['dx'],param['dy'])

dic = dict()
all_var2export = ['k','p']

for vars in all_var2export:
   exec('dic[vars] ='+vars)
    
np.save(runpath+'/analysis/spec_eke_b.npy',dic)
print('Everything stored.')

