from __future__ import print_function

path = '/network/aopp/cirrus/pred/kloewer/swm_back_ronew/'
outpath = '/network/home/aopp/kloewer/swm/paperplot/'

import numpy as np
import matplotlib.pyplot as plt
from cmocean import cm

plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'

runfolders = [0,3,8]

D = []
p = []

for r in runfolders:
    runpath = path+'run%04i' % r
    D.append(np.load(runpath+'/analysis/mean_e.npy').all()['em'])
    p.append(np.load(runpath+'/param.npy').all())

##
expo = 1/3.

for i in range(len(D)):
    D[i] = np.sign(D[i])*abs(D[i])**expo

## colormap

thermal2 = cm.thermal.from_list('thermal2',np.vstack(([0.5,0.5,0.5,1],cm.thermal(np.arange(1,256)))))
    

## plotting

levs = np.linspace(0,np.percentile(np.array(D),98),32)

s = 1e3

fig,axs = plt.subplots(1,len(D),sharex=True,sharey=True,figsize=(10,4.8))
fig.tight_layout(rect=[0.03,.15,1,0.9])
fig.subplots_adjust(wspace=0.03,hspace=0.03)

pos = axs[0].get_position()
pos2 = axs[-1].get_position()
cax = fig.add_axes([pos.x0,0.1,pos2.x1-pos.x0,0.03])

for i,(iD,ip) in enumerate(zip(D,p)):
    q = axs[i].contourf(ip['x_T']/s,ip['y_T']/s,iD,levs,extend='both',cmap=thermal2)
    
    
axs[0].set_title('LR + weak backscatter') 
axs[1].set_title('LR + moderate backscatter') 
axs[2].set_title('LR + strong backscatter') 

tiks = np.array([0,0.005,0.05,.2,.5,1,1.8,2.4])

cbar = fig.colorbar(q,cax=cax,orientation='horizontal',ticks=tiks**expo)
cbar.set_label(r'[m$^3$s$^{-2}$]')
cbar.set_ticklabels(tiks)

axs[0].set_xticks([0,1000,2000,3000])
axs[1].set_xticks([0,1000,2000,3000])
axs[2].set_xticks([0,1000,2000,3000])
axs[0].set_yticks([0,1000,2000,3000])

axs[0].set_title('a',loc='left',fontweight='bold')
axs[1].set_title('b',loc='left',fontweight='bold')
axs[2].set_title('c',loc='left',fontweight='bold')

axs[0].set_ylabel(r'$y$ [km]')
axs[1].set_xlabel(r'$x$ [km]')


plt.suptitle(r'Climatological mean subgrid-EKE $\overline{e}$',x=0.22)

plt.savefig(outpath+'plots/e_mean.png',dpi=300)
plt.close(fig)