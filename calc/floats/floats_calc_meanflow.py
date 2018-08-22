## LAGRANGIAN FLOATS
import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import RegularGridInterpolator as RGI

path = '/home/mkloewer/python/swm/'

# OPTIONS
runfolder = [14]
print('Calculating floats from run ' + str(runfolder))

## pad function for u,v
def uv_pad(u,v):
    """ Pads u,v with kinematic and lateral boundary conditions. """    
    bs = 0   # boundary speed (positive inwards)
    u = np.pad(u,((1,1),(1,1)),'constant',constant_values=((bs,-bs),(bs,-bs)))    # kinematic bc
    v = np.pad(v,((1,1),(1,1)),'constant',constant_values=((bs,-bs),(bs,-bs)))
        
    return u,v

runpath = path+'data/run%04i' % runfolder[0]

D = np.load(runpath+'/analysis/mean.npy').all()
param = np.load(runpath+'/param.npy').all()

u = D['um'].reshape((param['ny'],param['nx']-1))
v = D['vm'].reshape((param['ny']-1,param['nx']))
    
dt = 6*3600.

# renaming for convenience
x_u = param['x_u']
y_u = param['y_u']

x_v = param['x_v']
y_v = param['y_v']

# pad x,y to include the boundary conditions
# grid is not equidistant at the boundary however ...
x_u = np.hstack((0.,x_u,param['Lx']))
y_u = np.hstack((0.,y_u,param['Ly']))

y_v = np.hstack((0,y_v,param['Ly']))
x_v = np.hstack((0,x_v,param['Lx']))

print('Grid read.')

## OPTIONS
N = 10000     # number of floats
ntime = 1460    # number of timesteps to integrate forward: equals one year
R_chunk = 10

# preallocate
nbins = 254 
H = np.zeros((nbins,nbins))

#exclude up to 30km from boundary
Hrange = [[30e3,3810e3],[30e3,3810e3]]

# edges of inital seeding region
seedx = np.array([100,200])*1e3         
seedy = np.array([100,1920])*1e3
def u2mat(u):
    return u.reshape((param['ny'],param['nx']-1))

def v2mat(v):
    return v.reshape((param['ny']-1,param['nx']))

# padding with boundary conditions
u,v = uv_pad(u,v)
print('Padding boundary conditions: done.')

# Advect the floats
Xa = np.empty((ntime,N,R_chunk)).astype(np.float64)
Ya = np.empty_like(Xa)

# preallocate
X = np.empty((ntime,N)).astype(np.float64)
Y = np.empty_like(X).astype(np.float64)

for r in range(R_chunk):
    print(r)

    # set initial conditions
    X[0,:] = np.random.rand(N)*np.diff(seedx) + seedx[0]
    Y[0,:] = np.random.rand(N)*np.diff(seedy) + seedy[0]
    
    # interpolation function
    Iu = RGI((x_u,y_u),u.T,bounds_error=False,fill_value=0)
    Iv = RGI((x_v,y_v),v.T,bounds_error=False,fill_value=0)
    
    for i in range(ntime-1):
        # Following the ideas of the Crank-Nicolson method
        # However, as the RHS is discrete, Euler forward is used to estimate
        # the RHS of the next time step
    
        # old RHS
        dXdt1 = Iu((X[i,:],Y[i,:]))
        dYdt1 = Iv((X[i,:],Y[i,:]))
    
        # new RHS
        dXdt2 = Iu((X[i,:] + dt*dXdt1,Y[i,:] + dt*dYdt1))
        dYdt2 = Iv((X[i,:] + dt*dXdt1,Y[i,:] + dt*dYdt1))
        
        # average of old and new RHS
        X[i+1,:] = X[i,:] + 0.5*dt*(dXdt1+dXdt2)
        Y[i+1,:] = Y[i,:] + 0.5*dt*(dYdt1+dYdt2)
        
    Xa[:,:,r] = X
    Ya[:,:,r] = Y

print('Advecting floats done.')

# free memory
del u,v

# histogram do not count the grid cell at the boundary
Hi,xe,ye = np.histogram2d(Xa.flatten(),Ya.flatten(),bins=nbins,range=Hrange)
H += Hi
print('Histogram computed.')

# free memory
del Xa,Ya
    

print(H.sum())

# mid points
xm = xe[:-1] + (xe[1]-xe[0])/2.
ym = ye[:-1] + (ye[1]-ye[0])/2.

## STORING
dic = dict()
all_var2export = ['X','Y','seedx','seedy','H','xe','ye','xm','ym']

for vars in all_var2export:
    exec('dic[vars] ='+vars)
    
np.save(runpath+'/analysis/floats_meanflow.npy',dic)
print('Everything stored.')
