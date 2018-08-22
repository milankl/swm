## PRODUCE FINAL PLOTS
from __future__ import print_function

path1 = '/network/aopp/cirrus/pred/kloewer/swm_bf_cntrl/data/'
path2 = '/network/aopp/cirrus/pred/kloewer/swm_back_ronew/'
outpath = '/network/home/aopp/kloewer/swm/paperplot/'

import os; os.chdir(path2) # change working directory
import numpy as np
from scipy import sparse
from scipy.integrate import cumtrapz
import matplotlib.pyplot as plt
import time as tictoc
from netCDF4 import Dataset
import glob
from cmocean import cm

plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'

## OPTIONS
runfolders = [0,3,6]
print('Creating final plots from run ' + str(runfolders))

## read data
p = []
z = []

# import functions
funpath = '/network/home/aopp/kloewer/git/swm/'
exec(open(funpath+'swm_param.py').read())
exec(open(funpath+'swm_operators.py').read())

for i,r in enumerate(runfolders):
    
    path = path2 if i==1 else path1
    runpath = path+'run%04i' % r
    
    u = np.load(runpath+'/u_last.npy').flatten()
    v = np.load(runpath+'/v_last.npy').flatten()
    p.append(np.load(runpath+'/param.npy').all())
    
    # set param the current param dictionnary
    param = p[-1]
    param['output'] = 0
    
    # get operators
    set_grad_mat()
    set_interp_mat()
    set_coriolis()
    
    # rel vort calculation
    z.append(q2mat(Gvx.dot(v) - Guy.dot(u)))
    print((z[-1]**2).mean())
    expo = 0.5
    z[-1] = (np.sign(z[-1])*abs(z[-1]*1e5)**expo)

z_max = np.percentile(abs(z[-1]),98.5)
levs = np.linspace(-z_max,z_max,32)

##

fig,axs = plt.subplots(1,3,figsize=(10,4.8),sharex=True,sharey=True)
fig.tight_layout(rect=[0.03,.15,1,0.9])
fig.subplots_adjust(wspace=0.03,hspace=0.03)

pos = axs[0].get_position()
pos2 = axs[-1].get_position()
cax = fig.add_axes([pos.x0,0.1,pos2.x1-pos.x0,0.03])

axs[0].contourf(p[0]['x_q']/1e3,p[0]['y_q']/1e3,z[0],levs,cmap=cm.balance,extend='both')    
axs[1].contourf(p[1]['x_q']/1e3,p[1]['y_q']/1e3,z[1],levs,cmap=cm.balance,extend='both')    
q0 = axs[2].contourf(p[2]['x_q']/1e3,p[2]['y_q']/1e3,z[2],levs,cmap=cm.balance,extend='both')    

tik0 = np.array([-3,-1.9,-1,-0.4,-0.1,0,0.1,0.4,1,1.9,3])

cbar0 = fig.colorbar(q0,cax=cax,orientation='horizontal',ticks=np.sign(tik0)*abs(tik0)**expo)
cbar0.set_label(r'[10$^{-5}$ s$^{-1}$]')
cbar0.set_ticklabels(tik0)

plt.suptitle(r'Snapshots of relative vorticity',x=0.201)

axs[0].set_xlim(0,param['Lx']/1e3)
axs[0].set_ylim(0,param['Ly']/1e3)

axs[0].set_xticks([0,1000,2000,3000])
axs[0].set_yticks([0,1000,2000,3000])

axs[0].set_ylabel(r'$y$ [km]')
axs[1].set_xlabel(r'$x$ [km]')


axs[0].set_title(r'Low resolution, $\Delta x = $30km')
axs[1].set_title(r'LR + moderate backscatter')
axs[2].set_title(r'High resolution, $\Delta x = $7.5km')

axs[0].set_title('a',loc='left',fontweight='bold')
axs[1].set_title('b',loc='left',fontweight='bold')
axs[2].set_title('c',loc='left',fontweight='bold')

plt.savefig(outpath+'plots/relvort_snapshot.png',dpi=300)
plt.close(fig)
