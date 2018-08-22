## PLOT STREAMFUNCTION
from __future__ import print_function
#path = '/home/mkloewer/python/swm/'
path = '/network/aopp/cirrus/pred/kloewer/swm_bf_cntrl/data/'
path2 = '/home/kloewer/git/swm/'
outpath = '/network/home/aopp/kloewer/swm/paperplot/'

import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt
from cmocean import cm

plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'


# OPTIONS
runfolder = [0,6,8]
print('Compare mean plots from run ' + str(runfolder))

## read data

runpath1 = path+'run%04i' % runfolder[0]
D1 = np.load(runpath1+'/analysis/mean_pv.npy').all()
param1 = np.load(runpath1+'/param.npy').all()

runpath2 = path+'run%04i' % runfolder[1]
D2 = np.load(runpath2+'/analysis/mean_pv.npy').all()
param2 = np.load(runpath2+'/param.npy').all()

runpath3 = path+'run%04i' % runfolder[2]
D3 = np.load(runpath3+'/analysis/mean_pv.npy').all()
param3 = np.load(runpath3+'/param.npy').all()

def q2mat(q,param):
    return q.reshape((param['ny']+1,param['nx']+1))
    
def u2mat(u,param):
    return u.reshape((param['ny'],param['nx']-1))

def v2mat(v,param):
    return v.reshape((param['ny']-1,param['nx']))

def fq(param):
    y_q = param['y_q']
    Ly = param['Ly']

    # subtract the regions mid-y so that phi_0 corresponds to a central latitude
    yy_q = np.array([y_q - Ly/2.]*(param['nx']+1)).T

    # globally available coriolis parameters (only f_q is actually needed though)
    f_q = param['f_0'] + param['beta']*yy_q.flatten()
    
    return f_q

    
## Plotting
fig,axs = plt.subplots(2,3,sharex=True,sharey=True,figsize=(12,9))
fig.tight_layout(rect=[0.02,.15,1.,0.98])
fig.subplots_adjust(wspace=0.03,hspace=0.1)

pos = axs[0,0].get_position()
pos2 = axs[0,2].get_position()
cax = fig.add_axes([pos.x0,0.12,pos2.x1-pos.x0,0.02])
cax2 = fig.add_axes([pos.x0,0.05,pos2.x1-pos.x0,0.02])

s = 1e3 # scaling factor m -> km
n = 10

levs = np.linspace(6e-8,2.2e-7,17)
levs2 = np.linspace(-2e-8,2e-8,17)

# PV mean
q1 = axs[0,0].contourf(param1['x_q']/s,param1['y_q']/s,q2mat(D1['qm'],param1),levs,extend='both')
axs[0,1].contourf(param3['x_q']/s,param3['y_q']/s,q2mat(D3['qm'],param3),levs,extend='both')
axs[0,2].contourf(param2['x_q']/s,param2['y_q']/s,q2mat(D2['qm'],param2),levs,extend='both')

q2 = axs[1,0].contourf(param1['x_q']/s,param1['y_q']/s,q2mat(D1['qm']-fq(param1)/param1['H'],param1),levs2,extend='both',cmap=cm.balance)
axs[1,1].contourf(param3['x_q']/s,param3['y_q']/s,q2mat(D3['qm']-fq(param3)/param3['H'],param3),levs2,extend='both',cmap=cm.balance)
axs[1,2].contourf(param2['x_q']/s,param2['y_q']/s,q2mat(D2['qm']-fq(param2)/param2['H'],param2),levs2,extend='both',cmap=cm.balance)

# PV gradient
#axs[1].quiver(param3['x_q'][::n]/s,param3['y_q'][::n]/s,v2mat(D3['dqm_dx'],param3)[::n,::n],u2mat(D3['dqm_dy'],param3)[::n,::n])


fig.colorbar(q1,cax=cax,orientation='horizontal')
cbar = fig.colorbar(q2,cax=cax2,orientation='horizontal')
cbar.set_label(r'Potential Vorticity [m$^{-1}$s$^{-1}$]')
axs[0,0].set_xticks([0,1000,2000,3000])
axs[0,0].set_yticks([0,1000,2000,3000])

axs[0,0].set_xlim(0,param1['Lx']/s)
axs[0,0].set_ylim(0,param1['Ly']/s)
axs[0,0].set_ylabel(r'$y$ [km]')
axs[1,0].set_ylabel(r'$y$ [km]')
axs[1,0].set_xlabel(r'$x$ [km]')
axs[1,1].set_xlabel(r'$x$ [km]')#
axs[1,2].set_xlabel(r'$x$ [km]')

axs[0,0].set_title(r'PV, $\Delta x = $30km, no-slip')
axs[0,1].set_title(r'PV, $\Delta x = $15km, free-slip')
axs[0,2].set_title(r'PV, $\Delta x = $7.5km, no-slip')

axs[1,0].set_title(r'PV - f/H, $\Delta x = $30km, no-slip')
axs[1,1].set_title(r'PV - f/H, $\Delta x = $15km, free-slip')
axs[1,2].set_title(r'PV - f/H, $\Delta x = $7.5km, no-slip')

axs[0,0].set_title('a',fontweight='bold',loc='left')
axs[0,1].set_title('b',fontweight='bold',loc='left')
axs[0,2].set_title('c',fontweight='bold',loc='left')

axs[1,0].set_title('d',fontweight='bold',loc='left')
axs[1,1].set_title('e',fontweight='bold',loc='left')
axs[1,2].set_title('f',fontweight='bold',loc='left')

plt.savefig(outpath+'plots/pv_gradient.pdf')
plt.close(fig)
