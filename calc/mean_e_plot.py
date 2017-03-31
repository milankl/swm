from __future__ import print_function
path = '/home/mkloewer/python/swm/'
import os; os.chdir(path) # change working directory
import numpy as np
import matplotlib.pyplot as plt
from cmocean import cm

runfolders = [13,12,14]

D = []
p = []

for r in runfolders:
    runpath = path+'stoch/data/run%04i' % r
    D.append(np.load(runpath+'/analysis/mean_e.npy').all()['em'])
    p.append(np.load(runpath+'/param.npy').all())

print(np.percentile(np.array(D),95))
##
expo = (1/2.)

for i in range(len(D)):
    D[i] = np.sign(D[i])*abs(D[i]*1e4)**expo

## colormap

thermal2 = cm.thermal.from_list('thermal2',np.vstack(([0.5,0.5,0.5,1],cm.thermal(np.arange(1,256)))))
    

## plotting

levs = np.linspace(0,np.percentile(np.array(D),95),64)

fig,axs = plt.subplots(1,len(D),sharex=True,sharey=True,figsize=(10,5))
fig.tight_layout(rect=[0,.1,1,0.95])
fig.subplots_adjust(wspace=0.03,hspace=0.03)

pos = axs[0].get_position()
pos2 = axs[-1].get_position()
cax = fig.add_axes([pos.x0,0.1,pos2.x1-pos.x0,0.03])

for i,(iD,ip) in enumerate(zip(D,p)):
    q = axs[i].contourf(ip['x_T'],ip['y_T'],iD,levs,extend='both',cmap=thermal2)
    axs[i].set_title(r'$\overline{e}, n_{diss} = %.2f$' % ip['n_diss']) 

tiks = np.arange(0,15,2)**expo

cbar = fig.colorbar(q,cax=cax,orientation='horizontal',ticks=tiks)
cbar.set_label(r'$[cm^2/s^2]$')
cbar.set_ticklabels(np.round(tiks**2).astype(np.int))

axs[0].set_xticks([])
axs[0].set_yticks([])

plt.savefig(path+'compare/e_mean.png')
plt.close(fig)
#plt.show()