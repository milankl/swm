## COMPUTE ENSTROPHY SPECTRUM
from __future__ import print_function
path = '/home/mkloewer/python/swm/'
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
runfolder = [2,3]
print('Calculating EKE spectrogramms from run ' + str(runfolder))

##    
def ens_spec_avg(z,dx,dy):
    """ Computes a wavenumber-frequency plot for 3D (t,x,y) data via radial (k = sqrt(kx**2 + ky**2)) integration. TODO: correct normalisation, so that the integral in normal space corresponds to the integral in Fourier space.
    """
    
    nt,ny,nx = np.shape(z)
    kx = (1/(dx))*np.hstack((np.arange(0,(nx+1)/2.),np.arange(-nx/2.+1,0)))/float(nx)
    ky = (1/(dy))*np.hstack((np.arange(0,(ny+1)/2.),np.arange(-ny/2.+1,0)))/float(ny)

    kxx,kyy = np.meshgrid(kx,ky)
    # radial distance from kx,ky = 0
    kk = np.sqrt(kxx**2 + kyy**2) 

    if nx >= ny: #kill negative wavenumbers
        k  = kx[:int(nx/2)+1]
    else:
        k  = ky[:int(ny/2)+1]

    dk = k[1] - k[0]

    # create radial coordinates, associated with k[i]
    # nearest point interpolation to get points within the -.5,.5 annulus
    rcoords = []
    for i in range(len(k)):
        rcoords.append(np.where((kk>(k[i]-.5*dk))*(kk<=(k[i]+.5*dk))))

    # 2D FFT average
    
    pz = np.empty((nt,ny,nx))
    
    for i in range(nt):
        pz[i,:,:] = abs(np.fft.fft2(z[i,:,:]))**2
        if i % 100 == 0:
            print(i)
    
    pz_avg = .5*pz.mean(axis=0)

    # mulitply by dk to have the corresponding integral
    ens_spec = np.zeros(len(k))
    for i in range(len(k)):
        ens_spec[i] = np.sum(pz_avg[rcoords[i][0],rcoords[i][1]])*dk

    return k[1:],ens_spec[1:]  # eliminate zero wavenumber


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
tlen = len(time)
dt = time[1] - time[0]

## read param
global param
param = np.load(runpath+'/param.npy').all()
param['output'] = 0

set_grad_mat()
set_interp_mat()

#reshape u,v
u = u.reshape((tlen,param['Nu'])).T
v = v.reshape((tlen,param['Nv'])).T
z = (Gvx.dot(v) - Guy.dot(u)).T.reshape((tlen,param['ny']+1,param['nx']+1))
del u,v
##
k,p = ens_spec_avg(z,param['dx'],param['dy'])

##
dic = dict()
all_var2export = ['k','p']

for v in all_var2export:
    exec('dic[v] ='+v)
    
np.save(runpath+'/analysis/spec_ens.npy',dic)

print('Everything stored.')
