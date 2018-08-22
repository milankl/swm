## FLOATS PLOT
import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import RegularGridInterpolator as RGI
import matplotlib.pyplot as plt
from cmocean import cm

plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'

path = '/home/mkloewer/python/swm/'

# OPTIONS#
runfolder = [0,6,14]
print('Plotting floats from run ' + str(runfolder))

## LOAD DATA
runpath1 = path+'data/run%04i' % runfolder[0]
param1 = np.load(runpath1+'/param.npy').all()
D1 = np.load(runpath1+'/analysis/floats.npy').all()

runpath2 = path+'data/run%04i' % runfolder[1]
param2 = np.load(runpath2+'/param.npy').all()
D2 = np.load(runpath2+'/analysis/floats.npy').all()

runpath3 = path+'data/run%04i' % runfolder[2]
param3 = np.load(runpath3+'/param.npy').all()
D3 = np.load(runpath3+'/analysis/floats.npy').all()

expo = 0.5   # non-linear colorbar exponent

dx,dy = D1['xm'][1]-D1['xm'][0],D1['ym'][1]-D1['ym'][0]
D1['H'] = (D1['H']/dx/dy*1e6)**expo

dx,dy = D2['xm'][1]-D2['xm'][0],D2['ym'][1]-D2['ym'][0]
D2['H'] = (D2['H']/dx/dy*1e6)**expo

dx,dy = D1['xm'][1]-D1['xm'][0],D1['ym'][1]-D1['ym'][0]
D3['H'] = (D3['H']/dx/dy*1e6)**expo

tiks = np.array([0,1,5,10,20,30,40])#,50,60,70,80,100,120,140,160,180])
levs = np.linspace(0,np.percentile(D1['H'],99),64)

## Plotting
fig,axs = plt.subplots(2,3,sharex=True,sharey=True,figsize=(10,8))
plt.tight_layout(rect=[0,.09,1,0.98])
fig.subplots_adjust(wspace=0.03,hspace=0.03)

pos1 = axs[-1,0].get_position()
pos2 = axs[-1,-1].get_position()
cax = fig.add_axes([pos1.x0,0.07,pos2.x1-pos1.x0,0.03])

N = 30  # number of trajectories to plot

for iD,D in enumerate([D1,D3,D2]):
    for i in range(min(D['X'].shape[1],N)): # max num of floats to display
        z = np.random.randint(0,250,1)[0]
        axs[0,iD].plot(D['X'][:,40+i],D['Y'][:,40+i],color=cm.haline(z))

for i,D in enumerate([D1,D3,D2]):
    sx,sy = D['seedx'],D['seedy']   # unpack
    xm,ym,H = D['xm'],D['ym'],D['H']
    axs[1,i].add_patch(plt.Rectangle((sx[0],sy[0]),sx[1]-sx[0],sy[1]-sy[0],fc='none',alpha=1,ec='k',lw=1.5))
    q = axs[1,i].contourf(xm,ym,H.T,levs,extend='max',cmap=cm.thermal)
    axs[1,i].contour(xm,ym,H.T,[0],colors='w')

cbar = plt.colorbar(q,cax=cax,orientation='horizontal',ticks=tiks**expo)
cbar.set_ticklabels(tiks)
cbar.set_label('[1/km$^2$]')

axs[0,0].set_ylabel('Example float trajectories')
axs[1,0].set_ylabel('Accumulated float density')
axs[0,0].set_xlim(0,param1['Lx'])
axs[0,0].set_ylim(0,param1['Ly'])
axs[0,0].set_xticks([])
axs[0,0].set_yticks([])

axs[0,0].set_title(r'Low resolution, $\Delta x = $30km')
axs[0,2].set_title(r'High resolution, $\Delta x = $7.5km')
axs[0,1].set_title(r'LR +  moderate backscatter')

abc = 'abcdef'
abcc = 'kkkwww'
abci = 0
for axcol in axs:
    for ax in axcol:
        plt.text(0.93,0.05,abc[abci],transform=ax.transAxes,fontweight='bold',color=abcc[abci])
        abci += 1


axs[-1,0].set_xlabel(r'$x$')
axs[-1,1].set_xlabel(r'$x$')
axs[-1,2].set_xlabel(r'$x$')

axs[0,2].set_ylabel(r'$y$')
axs[0,2].yaxis.set_label_position('right')

axs[1,2].set_ylabel(r'$y$')
axs[1,2].yaxis.set_label_position('right')

plt.savefig(path+'compare/floats_3runs_bf.png',dpi=150)
plt.close(fig)