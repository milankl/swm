## LAGRANGIAN FLOATS
import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import RegularGridInterpolator as RGI

path = '/home/mkloewer/python/swm/'

# OPTIONS
runfolder = [2,3,4,5,6]
print('Calculating floats from run ' + str(runfolder))

## pad function for u,v
def uv_pad(u,v):
    """ Pads u,v with kinematic and lateral boundary conditions. """    
    bs = 0   # boundary speed (positive inwards)
    u = np.pad(u,((0,0),(1,1),(1,1)),'constant',constant_values=((0,0),(bs,-bs),(bs,-bs)))    # kinematic bc
    v = np.pad(v,((0,0),(1,1),(1,1)),'constant',constant_values=((0,0),(bs,-bs),(bs,-bs)))
        
    return u,v

## Reading data in blocks to save memory.
def read_data(index_start):
    """ Based on the start index decide which datasets to load and do so."""
    index_end = index_start + 2*ntime
    
    file_length0 = np.hstack((0,file_length))
    file_length_each = np.hstack((file_length[0],np.diff(file_length)))
    
    start_crit = np.logical_and(index_start <= file_length,index_start > file_length0[:-1])
    end_crit = np.logical_and(index_end <= file_length,index_end > file_length0[:-1])
    
    data_to_read = np.logical_or(start_crit,end_crit)
    
    d = 0
    for r,rn in zip(np.array(runfolder)[data_to_read],np.arange(len(runfolder))[data_to_read]):
        
        # r is the run id
        # rn is the relative number of the run, starting with 0 for the first in runfolder
        
        runpath = path+'data/run%04i' % r
        
        if d == 0: # in case the first data set is read
            
            idx0 = index_start - file_length0[rn]
            idx1 = min(index_end - file_length0[rn],file_length_each[rn]+1)

            ncu = Dataset(runpath+'/u.nc')
            u = ncu.variables['u'][idx0:idx1,:,:]
            ncu.close()
    
            ncv = Dataset(runpath+'/v.nc')
            v = ncv.variables['v'][idx0:idx1,:,:]
            ncv.close()
    
        else:   # concatenate second data set
            print('reading second file')
            idx0 = 1    # skip the first as this is a duplicate from the last step of previous file
            idx1 = index_end - file_length0[rn] + 1
        
            ncu = Dataset(runpath+'/u.nc')
            ncv = Dataset(runpath+'/v.nc')
            
            u = np.concatenate((u,ncu.variables['u'][idx0:idx1,:,:]))
            v = np.concatenate((v,ncv.variables['v'][idx0:idx1,:,:]))
        
        d = 1   # set to 1 to concatenate any further data set
    
    return u,v

## read grid
file_length = np.empty(len(runfolder)).astype(np.int)   # list to store the number of time steps in each file

for i,r in enumerate(runfolder):
    runpath = path+'data/run%04i' % r
    ncu = Dataset(runpath+'/u.nc')
    
    if i == 0: time = ncu.variables['t'][:]
    else: time = np.hstack((time,ncu.variables['t'][:]))
    
    file_length[i] = len(time)
    
dt = time[1]-time[0]

# read spatial grid
param = np.load(runpath+'/param.npy').all()

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
N = 100         # number of floats
ntime = 1460    # number of timesteps to integrate forward: equals one year
R = 1000        # number of repetitions (from random start dates)

# preallocate
nbins = 254 
H = np.zeros((nbins,nbins))

#exclude up to 30km from boundary
Hrange = [[30e3,3810e3],[30e3,3810e3]]

# edges of inital seeding region
seedx = np.array([100,200])*1e3         
seedy = np.array([100,1920])*1e3

all_startidxs = np.sort(np.random.randint(0,len(time)-ntime-2,R))

ichunk = 0
while all_startidxs.shape[0]:   # delete from all_startidxs in every loop until empty
    ichunk += 1

    # indices for current chunk
    startidxs = all_startidxs[all_startidxs < (all_startidxs[0]+ntime)]
    
    # number of repetitions within chunk
    R_chunk = len(startidxs)    
    
    # remove these dates from all_startidxs
    all_startidxs = all_startidxs[len(startidxs):]

    u,v = read_data(startidxs[0])
    print('Data chunk %i read.' % ichunk)
    
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
        # relative index within dataset
        i_rel = startidxs[r] - startidxs[0]
    
        # set initial conditions
        X[0,:] = np.random.rand(N)*np.diff(seedx) + seedx[0]
        Y[0,:] = np.random.rand(N)*np.diff(seedy) + seedy[0]
        
        # old time step
        Iu1 = RGI((x_u,y_u),u[i_rel,:,:].T,bounds_error=False,fill_value=-1)
        Iv1 = RGI((x_v,y_v),v[i_rel,:,:].T,bounds_error=False,fill_value=-1)
        
        for i in range(ntime-1):
            # Following the ideas of the Crank-Nicolson method
            # However, as the RHS is discrete, Euler forward is used to estimate
            # the RHS of the next time step
            
            # next time step
            # fill_value means once they are out of the boundary floats stop moving
            Iu2 = RGI((x_u,y_u),u[i+i_rel+1,:,:].T,bounds_error=False,fill_value=0)
            Iv2 = RGI((x_v,y_v),v[i+i_rel+1,:,:].T,bounds_error=False,fill_value=0)
        
            # old RHS
            dXdt1 = Iu1((X[i,:],Y[i,:]))
            dYdt1 = Iv1((X[i,:],Y[i,:]))
        
            # new RHS based on new time step
            dXdt2 = Iu2((X[i,:] + dt*dXdt1,Y[i,:] + dt*dYdt1))
            dYdt2 = Iv2((X[i,:] + dt*dXdt1,Y[i,:] + dt*dYdt1))
            
            # average of old and new RHS
            X[i+1,:] = X[i,:] + 0.5*dt*(dXdt1+dXdt2)
            Y[i+1,:] = Y[i,:] + 0.5*dt*(dYdt1+dYdt2)
            
            # rename to avoid re-setting up of RGI functions
            Iu1 = Iu2   # next -> old
            Iv1 = Iv2
        
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
    
np.save(runpath+'/analysis/floats.npy',dic)
print('Everything stored.')
