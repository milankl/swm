## PLOT ENERGY AND ENSTROPY TIME SERIES
from __future__ import print_function

path1 = '/network/aopp/cirrus/pred/kloewer/swm_back_ronew/'
path2 = '/network/aopp/cirrus/pred/kloewer/swm_bf_cntrl/data/'
outpath = '/network/home/aopp/kloewer/swm/paperplot/'

import os; os.chdir(path2) # change working directory
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec

plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'

# OPTIONS
runfolder = [0,6,0,3,8]
print('Compare mean plots from run ' + str(runfolder))

## read data

runpath1 = path2+'run%04i' % runfolder[0]
D1 = np.load(runpath1+'/analysis/mean_timeseries.npy').all()

runpath2 = path2+'run%04i' % runfolder[1]
D2 = np.load(runpath2+'/analysis/mean_timeseries.npy').all()

runpath3 = path1+'run%04i' % runfolder[2]
D3 = np.load(runpath3+'/analysis/mean_timeseries.npy').all()

runpath4 = path1+'run%04i' % runfolder[3]
D4 = np.load(runpath4+'/analysis/mean_timeseries.npy').all()

runpath5 = path1+'run%04i' % runfolder[4]
D5 = np.load(runpath5+'/analysis/mean_timeseries.npy').all()

for D in [D1,D3,D4,D5]:  # from seconds to years
    D['t'] = D['t'] / 3600. / 24. / 365.

# create new t2 time axis based on the knowledge that dt = 86260 seconds
D2['t'] = (np.arange(len(D2['t']))*86260) / 3600. / 24. / 365.

## 

def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0)) 
    return (cumsum[N:] - cumsum[:-N]) / N 

## PLOTTING   

#fig,((ax1,ax3),(ax2,ax4)) = plt.subplots(2,2,sharex=True,sharey='row',figsize=(8,6))

#fig = plt.figure(figsize=(8,6))parameteri
#gs = gridspec.GridSpec(2, 2,width_ratios=[20,1])

# fig = plt.figure(figsize=(8,3))
# gs = gridspec.GridSpec(1, 2,width_ratios=[20,1])

fig,ax1 = plt.subplots(1,1,figsize=(8,3))

# ax1 = plt.subplot(gs[0])
#ax2 = plt.subplot(gs[1,0])

# ax3 = plt.subplot(gs[1])
#ax4 = plt.subplot(gs[1,1])

fig.tight_layout(rect=[0.02,0.03,1.01,0.95])
# fig.subplots_adjust(wspace=0.03,hspace=0.15)

#ax1.set_xticklabels([])

# ax3.set_xticks([])
# ax4.set_xticks([])

# ax3.set_yticklabels([])
# ax4.set_yticklabels([])

ax1.set_xlim(0,30)
# ax2.set_xlim(0,30)
ax1.set_ylim(0,170)
# ax3.set_ylim(0,170)

# ax2.set_ylim(0,10)
# ax4.set_ylim(0,10)

# ax3.set_xlim(0,1)
# ax4.set_xlim(0,1)

KEmm = [0,]*5
PEmm = [0,]*5
for i,D in enumerate([D1,D2,D3,D4,D5]):
    KEmm[i] = D['KEm'][5*365:].mean()/1e3
    PEmm[i] = D['PEm'][5*365:].mean()/1e3

# ax3.plot([0,1],np.ones(2)*KEmm[1],'C2',lw=3)
# ax3.plot([0,1],np.ones(2)*KEmm[0],'C0',lw=3)
# ax3.plot([0,1],np.ones(2)*KEmm[2],'C1',lw=1)
# ax3.plot([0,1],np.ones(2)*KEmm[3],'C3',lw=1)
# ax3.plot([0,1],np.ones(2)*KEmm[4],'C5',lw=1)

# ax3.text(0.2,KEmm[1],"%i" % KEmm[1])
# ax3.text(0.2,KEmm[0],"%i" % KEmm[0])
# ax3.text(0.2,KEmm[2],"%i" % KEmm[2])
# ax3.text(0.2,KEmm[3],"%i" % KEmm[3])
# ax3.text(0.2,KEmm[4],"%i" % KEmm[4])

# ax4.plot([0,1],np.ones(2)*PEmm[1],'C2',lw=3)
# ax4.plot([0,1],np.ones(2)*PEmm[0],'C0',lw=3)
# ax4.plot([0,1],np.ones(2)*PEmm[2],'C1',lw=1)
# ax4.plot([0,1],np.ones(2)*PEmm[3],'C3',lw=1)
# ax4.plot([0,1],np.ones(2)*PEmm[4],'C5',lw=1)

# ax4.text(0.2,PEmm[1],"%.1f" % PEmm[1])
# ax4.text(0.2,PEmm[0],"%.1f" % PEmm[0])
# ax4.text(0.2,PEmm[2],"%.1f" % PEmm[2])
# ax4.text(0.2,PEmm[3],"%.1f" % PEmm[3])
# ax4.text(0.2,PEmm[4],"%.1f" % PEmm[4])

#
ax1.plot(D2['t'],D2['KEm']/1e3,'C2',label=r'High resolution, $\Delta x = $7.5km',lw=1.5)
ax1.plot(D1['t'],D1['KEm']/1e3,'C0',label=r'Low resolution, $\Delta x = $30km',lw=1.5)
ax1.plot(D3['t'],D3['KEm']/1e3,'C1',label=r'LR + weak backscatter',lw=1.5)
ax1.plot(D4['t'],D4['KEm']/1e3,'C3',label=r'LR + moderate backscatter',lw=1.5)
ax1.plot(D5['t'],D5['KEm']/1e3,'C5',label=r'LR + strong backscatter',lw=1.5)

# monthly mean for potential energy
m = 1
for D in [D1,D2,D3,D4,D5]:
    #D['PEm'] = D['PEm'][:-(len(D['PEm']) % m)].reshape((-1,m)).mean(axis=1)
    #D['t'] = D['t'][:-(len(D['t']) % m)].reshape((-1,m)).mean(axis=1)
    
    D['PEm'] = running_mean(D['PEm'],m)
    D['t'] = running_mean(D['t'],m)

# ax2.plot(D2['t'],D2['PEm']/1e3,'C2',lw=2)
# ax2.plot(D1['t'],D1['PEm']/1e3,'C0',lw=2)
# ax2.plot(D3['t'],D3['PEm']/1e3,'C1',lw=1)
# ax2.plot(D4['t'],D4['PEm']/1e3,'C3',lw=1)
# ax2.plot(D5['t'],D5['PEm']/1e3,'C5',lw=1,zorder=-1)

ax1.add_patch(patches.Rectangle((0,0),5,ax1.get_ylim()[1],color='0.7'))
# ax2.add_patch(patches.Rectangle((0,0),5,ax2.get_ylim()[1],color='grey',alpha=.3))
#invisible bar for legend
ax1.bar(0,0,color='0.7',label='spin-up')

ax1.set_title('Kinetic energy',loc='left')
#ax1.set_title('a',loc='left',fontweight='bold')
ax1.set_ylabel(r'kJ m$^{-2}$')
ax1.set_ylim(0,150)

# ax2.set_title('Potential energy PE')
# ax2.set_title('b',loc='left',fontweight='bold')
# ax2.set_ylabel(r'kJm$^{-2}$')
ax1.set_xlabel('time [years]')
ax1.legend(loc=4,fontsize=8,ncol=3)

#ax1.grid()

# ax3.set_title('c',fontweight='bold')
# ax4.set_title('d',fontweight='bold')

# ax3.set_ylabel('mean(KE)')
# ax3.yaxis.set_label_position('right')

# ax4.set_ylabel('mean(PE)')
# ax4.yaxis.set_label_position('right')


plt.savefig(outpath + 'plots/energy_timeseries_ke.eps')
plt.close(fig)
