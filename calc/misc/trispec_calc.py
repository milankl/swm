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
runfolder = [3]
print('Calculating 3D spectrogramms from run ' + str(runfolder))

##
def trispec(a,dt,dx,dy):
    """ Computes a wavenumber-frequency plot for 3D (t,x,y) data via radial (k = sqrt(kx**2 + ky**2)) integration. TODO: correct normalisation, so that the integral in normal space corresponds to the integral in Fourier space.
    """
    
    nt,ny,nx = np.shape(a)
    kx = (1/(dx))*np.hstack((np.arange(0,(nx+1)/2.),np.arange(-nx/2.+1,0)))/float(nx)
    ky = (1/(dy))*np.hstack((np.arange(0,(ny+1)/2.),np.arange(-ny/2.+1,0)))/float(ny)
    f = (1/(dt))*np.hstack((np.arange(0,(nt+1)/2.),np.arange(-nt/2.+1,0)))/float(nt)

    kxx,kyy = np.meshgrid(kx,ky)
    # radial distance from kx,ky = 0
    kk = np.sqrt(kxx**2 + kyy**2) 

    if nx >= ny: #kill negative wavenumbers
        k  = kx[:int(nx/2)+1]
    else:
        k  = ky[:int(ny/2)+1]

    f = f[:int(nt/2)+1] #kill negative frequencies
    dk = k[1] - k[0]

    # create radial coordinates, associated with k[i]
    # nearest point interpolation to get points within the -.5,.5 annulus
    rcoords = []
    for i in range(len(k)):
        rcoords.append(np.where((kk>(k[i]-.5*dk))*(kk<=(k[i]+.5*dk))))

    # 3D FFT
    p = np.fft.fftn(a)
    p = np.real(p * np.conj(p))

    # mulitply by dk to have the corresponding integral
    spec = np.zeros((len(f),len(k)))
    for i in range(len(k)):
        spec[:,i] = np.sum(p[:int(nt/2)+1,rcoords[i][0],rcoords[i][1]],axis=1)*dk

    return k[1:],f[1:],spec[1:,1:]  # eliminate zero frequency and zero wavenumber

## read data
for r,i in zip(runfolder,range(len(runfolder))):
    runpath = path+'data/run%04i' % r
    
    if i == 0:
        #u = np.load(runpath+'/u_sub.npy')
        #v = np.load(runpath+'/v_sub.npy')
        h = np.load(runpath+'/h_sub.npy')        
        time = np.load(runpath+'/t_sub.npy')
        print('run %i read.' % r)

    else:
        #u = np.concatenate((u,np.load(runpath+'/u_sub.npy')))
        #v = np.concatenate((v,np.load(runpath+'/v_sub.npy')))
        h = np.concatenate((h,np.load(runpath+'/h_sub.npy')))
        time = np.hstack((time,np.load(runpath+'/t_sub.npy')))
        print('run %i read.' % r)

t = time / 3600. / 24.  # in days
tlen = len(time)
dt = time[1] - time[0]

## read param
global param
param = np.load(runpath+'/param.npy').all()

##
k,f,p = trispec(h,dt,param['dx'],param['dy'])

##
dic = dict()
all_var2export = ['f','k','p']

for v in all_var2export:
    exec('dic[v] ='+v)
    
np.save(runpath+'/analysis/trispec_h.npy',dic)
print('Everything stored.')