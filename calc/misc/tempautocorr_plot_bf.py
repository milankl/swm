## TEMPAUTOCORR PLOT
from __future__ import print_function
path = '/home/mkloewer/python/swm/'
import os; os.chdir(path) # change working directory
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumtrapz

## read data

runpath1 = path+'data/run%04i' % 0
D1 = np.load(runpath1+'/analysis/acfs.npy').all()
param1 = np.load(runpath1+'/param.npy').all()

runpath2 = path+'data/run%04i' % 6
D2 = np.load(runpath2+'/analysis/acfs.npy').all()
param2 = np.load(runpath2+'/param.npy').all()

runpath3 = path+'data/run%04i' % 10
D3 = np.load(runpath3+'/analysis/acfs.npy').all()
param3 = np.load(runpath3+'/param.npy').all()

runpath4 = path+'data/run%04i' % 14
D4 = np.load(runpath4+'/analysis/acfs.npy').all()
param4 = np.load(runpath4+'/param.npy').all()

runpath5 = path+'data/run%04i' % 15
D5 = np.load(runpath5+'/analysis/acfs.npy').all()
param5 = np.load(runpath5+'/param.npy').all()

# read without bottom friction data

# runpath1 = path+'data/newold/run%04i' % 3
# D1 = np.load(runpath1+'/analysis/acfs.npy').all()
# param1 = np.load(runpath1+'/param.npy').all()
# 
# runpath2 = path+'data/newold/run%04i' % 10
# D2 = np.load(runpath2+'/analysis/acfs.npy').all()
# param2 = np.load(runpath2+'/param.npy').all()
# 
# runpath3 = path+'stoch/data/run%04i' % 13
# D3 = np.load(runpath3+'/analysis/acfs.npy').all()
# param3 = np.load(runpath3+'/param.npy').all()
# 
# runpath4 = path+'stoch/data/run%04i' % 12
# D4 = np.load(runpath4+'/analysis/acfs.npy').all()
# param4 = np.load(runpath4+'/param.npy').all()
# 
# runpath5 = path+'stoch/data/run%04i' % 14
# D5 = np.load(runpath5+'/analysis/acfs.npy').all()
# param5 = np.load(runpath5+'/param.npy').all()

## Plotting

fig,axs = plt.subplots(2,3,sharex=True,sharey=True,figsize=(9,6))
plt.tight_layout(rect=[0.05,0.02,1,0.96])
fig.subplots_adjust(wspace=0.05,hspace=0.26)

for i in range(3):
    for j in range(2):
        axs[j,i].plot(D1['time'],D1['acfs'][:,i,j],'C0',label=r'Low resolution, $\Delta x = $30km',lw=2)
        axs[j,i].plot(D2['time'],D2['acfs'][:,i,j],'C2',label=r'High resolution, $\Delta x = $7.5km',lw=2)
        axs[j,i].plot(D3['time'],D3['acfs'][:,i,j],'C3',label=r'LR + weak backscatter',lw=1,ls='--')
        axs[j,i].plot(D4['time'],D4['acfs'][:,i,j],'C1',label=r'LR + moderate backscatter',lw=1,ls='--')
        axs[j,i].plot(D5['time'],D5['acfs'][:,i,j],'C5',label=r'LR + strong backscatter',lw=1,ls='--')
        axs[j,i].plot([0,50],[0,0],'C7',alpha=.5)
    
axs[0,0].set_xlim(0,49)
axs[0,0].set_ylim(-0.9,1)

axs[1,1].legend(bbox_to_anchor=(-1, 1.02, 3., .102), loc=3, fontsize=9, ncol=3, mode="expand", borderaxespad=0.)

axs[0,0].set_title(r'Zonal velocity $u$',loc='left')
axs[0,1].set_title(r'Meridional velocity $v$',loc='left')
axs[0,2].set_title(r'Surface displacement $\eta$',loc='left')

axs[1,0].set_xlabel('lag [days]')
axs[1,1].set_xlabel('lag [days]')
axs[1,2].set_xlabel('lag [days]')

axs[0,0].set_ylabel('Autocorrelation \n Point A')
axs[1,0].set_ylabel('Autocorrelation \n Point B')

axs[0,0].set_yticks([-0.5,0,0.5,1])

abc = 'abcdef'
abci = 0
for axcol in axs:
    for ax in axcol:
        plt.text(0.93,0.93,abc[abci],transform=ax.transAxes,fontweight='bold')
        abci += 1

plt.savefig(path+'compare/autocorrelation_bf.pdf')
plt.close(fig)
#plt.show()