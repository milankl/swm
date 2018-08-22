import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'

param = dict()

param['nx'] = 128               # number of grid points in x-direction
param['ny'] = 128               # number of grid points in y-direction

param['Lx'] = 3840e3            # basin width L [meters]
param['Ly'] = 3840e3            # north-south basin extent [meters]

param['dx'] = param['Lx'] / param['nx'] # grid spacing in x-direction
param['dy'] = param['Ly'] / param['ny'] # grid spacing in y-direction

param['x_T'] = np.arange(param['dx']/2.,param['Lx'],param['dx'])
param['y_T'] = np.arange(param['dy']/2.,param['Ly'],param['dy'])

# grid vectors for u-points
param['x_u'] = param['x_T'][:-1] + param['dx']/2.
param['y_u'] = param['y_T']

Lx,Ly = param['Lx'],param['Ly']     # for convenience
param['rho'] = 1e3                  # density

xx_u,yy_u = np.meshgrid(param['x_u'],param['y_u'])

param['Fx0'] = 0.12      # was 0.12
Fx = param['Fx0']*(np.cos(2*np.pi*(yy_u-Ly/2)/Ly) + 2*np.sin(np.pi*(yy_u - Ly/2)/Ly))


fig = plt.figure(figsize=(6, 4.5)) 
gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1]) 
ax0 = plt.subplot(gs[0])
ax1 = plt.subplot(gs[1])

s = 6
ax0.quiver(xx_u[2::s,::s]/1e3,yy_u[2::s,::s]/1e3,Fx[2::s,::s],np.zeros_like(Fx)[2::s,::s], scale=6, headwidth=5)

ax0.set_title('a',loc='left',fontweight='bold')
ax0.set_title(r'$(F_x,F_y)$') 
ax0.set_xlabel(r'$x$ [km]')
ax0.set_ylabel(r'$y$ [km]')

ax0.set_xlim(0,Lx/1e3)
ax0.set_ylim(0,Ly/1e3)

ax0.scatter(300,1920,40,'C0',label='Point A')
ax0.scatter(2880,1920,40,'C1',label='Point B')

ax0.legend(loc=1)

ax1.plot(Fx[:,0],param['y_u'],lw=1.5)
ax1.plot(np.zeros_like(param['y_u']),param['y_u'],lw=0.5,alpha=0.5,color='k')
ax1.set_title('b',loc='left',fontweight='bold')
ax1.set_title(r'$\gamma(y)$')
ax1.set_ylim(0,Ly)
ax1.set_yticks([])
ax1.set_xticks([-0.2,0,0.2])
ax1.set_xticklabels([-0.2,0,0.2])
ax1.set_xlabel('[Pa]')




plt.tight_layout()
plt.savefig('/home/mkloewer/Dropbox/thesis/Chapter1/Chapter1Figs/forcing.pdf')