## FLOATS PLOT
import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import RegularGridInterpolator as RGI
import matplotlib.pyplot as plt
from cmocean import cm

path = '/home/mkloewer/python/swm/'

# OPTIONS#
runfolder = [3,10,7]
print('Plotting floats from run ' + str(runfolder))

## LOAD DATA
runpath1 = path+'data/run%04i' % runfolder[0]
param1 = np.load(runpath1+'/param.npy').all()
D1 = np.load(runpath1+'/analysis/floats.npy').all()

runpath2 = path+'data/run%04i' % runfolder[1]
param2 = np.load(runpath2+'/param.npy').all()
D2 = np.load(runpath2+'/analysis/floats.npy').all()

runpath3 = path+'data/run%04i' % runfolder[1]
param3 = np.load(runpath3+'/param.npy').all()
D3 = np.load(runpath3+'/analysis/floats.npy').all()

## Plotting
fig,axs = plt.subplots(2,3,sharex=True,sharey=True,figsize=(10,8))
plt.tight_layout(rect=[0,.07,1,0.98])
fig.subplots_adjust(wspace=0.03,hspace=0.03)

pos1 = axs[-1,0].get_position()
pos2 = axs[-1,-1].get_position()
cax = fig.add_axes([pos1.x0,0.07,pos2.x1-pos1.x0,0.03])

N = 30  # number of trajectories to plot

for iD,D in enumerate([D1,D2,D3]):
    for i in range(min(D['X'].shape[1],N)): # max num of floats to display
        z = np.random.randint(0,250,1)[0]
        axs[0,iD].plot(D['X'][:,i],D['Y'][:,i],color=cm.haline(z))

for i,D in enumerate([D1,D2,D3]):
    sx,sy = D['seedx'],D['seedy']   # unpack
    xm,ym,H = D['xm'],D['ym'],D['H']
    dx,dy = xm[1]-xm[0],ym[1]-ym[0]
    print(H.sum()/1460/100/1000)
    H = H/dx/dy*1e6 # float density per square kilometer
    axs[1,i].add_patch(plt.Rectangle((sx[0],sy[0]),sx[1]-sx[0],sy[1]-sy[0],fc='k',alpha=0,ec='k',lw=3))
    axs[1,i].add_patch(plt.Rectangle((sx[0],sy[0]),sx[1]-sx[0],sy[1]-sy[0],fc='none',alpha=1,ec='k',lw=1.5))
    q = axs[1,i].contourf(xm,ym,H.T,np.linspace(0,40,64),extend='max',cmap=cm.thermal)
    axs[1,i].contour(xm,ym,H.T,[0],colors='w')

cbar = plt.colorbar(q,cax=cax,orientation='horizontal',ticks=[0,10,20,30,40])
cbar.set_label('[1/km$^2$]')

axs[0,0].set_ylabel('Float trajectories (N=%i)' % N)
axs[1,0].set_ylabel('Accumulated float density (N=1e5)')
axs[0,0].set_xlim(0,param1['Lx'])
axs[0,0].set_ylim(0,param1['Ly'])
axs[0,0].set_xticks([])
axs[0,0].set_yticks([])

axs[0,0].set_title(r'Low resolution, $\Delta x = 30$km')
axs[0,1].set_title(r'High resolution, $\Delta x = 7.5$km')
axs[0,2].set_title(r'Low resolution + backscatter')

plt.savefig(path+'compare/floats_3runs.png')
plt.close(fig)