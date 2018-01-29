## COMPUTE EKE SPECTRUM
from __future__ import print_function

# path
import os
#path = os.path.dirname(os.getcwd()) + '/'   # on level above
path = '/network/aopp/cirrus/pred/kloewer/swm_back_ronew/'
os.chdir(path)                              # change working directory

import numpy as np
from scipy import sparse

# OPTIONS
# several entries in the list concatenates the runs and stores the result in the last folder
runfolder = [5]     
print('Calculating eke spectrum from run ' + str(runfolder))

##
def eke_spec_avg(u,v,dx,dy):
    """ Computes a wavenumber-frequency plot for 3D (t,x,y) data via radial (k = sqrt(kx**2 + ky**2)) integration. TODO: correct normalisation, so that the integral in normal space corresponds to the integral in Fourier space.
    """
    
    nt,ny,nx = np.shape(u)
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
    rcoords = []
    for i in range(len(k)):
        rcoords.append(np.where(kk<k[i]))

    # mulitply by dk to have the corresponding integral
    eke_spec = np.zeros(len(k))
    for i in range(len(k)):
        eke_spec[i] = np.sum(p_eke_avg[rcoords[i][0],rcoords[i][1]])
    
    eke_spec = np.diff(eke_spec) / dk
    k = (k[:-1] + k[1:])/2.

    return k,eke_spec

## read data
for r,i in zip(runfolder,range(len(runfolder))):
    runpath = path+'run%04i' % r
    
    if i == 0:
        skip = 5*365
        u = np.load(runpath+'/u_sub.npy')[skip:,...]
        v = np.load(runpath+'/v_sub.npy')[skip:,...]
        t = np.load(runpath+'/t_sub.npy')[skip:,...]
        print('run %i read.' % r)

    else:
        u = np.concatenate((u,np.load(runpath+'/u_sub.npy')))
        v = np.concatenate((v,np.load(runpath+'/v_sub.npy')))
        t = np.hstack((t,np.load(runpath+'/t_sub.npy')))
        print('run %i read.' % r)

## read param
global param
param = np.load(runpath+'/param.npy').all()

tlen = len(t)
dt = t[1] - t[0]

## create ouputfolder
try:
    os.mkdir(runpath+'/analysis')
except:
   pass

## include boundary conditions
u = np.pad(u,((0,0),(0,0),(1,1)),'constant')    # kinematic
v = np.pad(v,((0,0),(1,1),(0,0)),'constant')

u = np.pad(u,((0,0),(0,1),(0,0)),'edge')    # first: free-slip
v = np.pad(v,((0,0),(0,0),(0,1)),'edge')

u[:,-1,:] = (1-param['lbc'])*u[:,-1,:]  # now: adapt the actual boundary condition
v[:,:,-1] = (1-param['lbc'])*v[:,:,-1]

## actual calculation
k,p = eke_spec_avg(u,v,param['dx'],param['dy'])

## storing
dic = dict()
all_var2export = ['k','p']

for vars in all_var2export:
   exec('dic[vars] ='+vars)
    
np.save(runpath+'/analysis/spec_eke.npy',dic)
print('Everything stored.')

