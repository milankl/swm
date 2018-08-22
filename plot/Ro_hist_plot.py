## HISTOGRAM PLOTTING FOR REYNOLDS AND ROSSBY NUMBERS
from __future__ import print_function

path = '/network/home/aopp/kloewer/swm/paperplot/'
path2 = '/network/aopp/cirrus/pred/kloewer/swm_bf_cntrl/'
path3 = '/network/aopp/cirrus/pred/kloewer/swm_back_ronew/'

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'

# OPTIONS
runfolder = [0,6,0,2,4]
print('Creating Ro histogram plot from run '+str(runfolder))

## read data
runpath = path2+'data/run%04i' % runfolder[0]
D1 = np.load(runpath+'/analysis/Ro_defo_hist.npy').all()

runpath = path2+'data/run%04i' % runfolder[1]
D2 = np.load(runpath+'/analysis/Ro_defo_hist.npy').all()

runpath = path3+'run%04i' % runfolder[2]
D3 = np.load(runpath+'/analysis/Ro_defo_hist.npy').all()
print(np.load(runpath+'/param.npy').all()['n_diss'])

runpath = path3+'run%04i' % runfolder[3]
D4 = np.load(runpath+'/analysis/Ro_defo_hist.npy').all()
print(np.load(runpath+'/param.npy').all()['n_diss'])

runpath = path3+'run%04i' % runfolder[4]
D5 = np.load(runpath+'/analysis/Ro_defo_hist.npy').all()
print(np.load(runpath+'/param.npy').all()['n_diss'])

## EXTEND HISTOGRAM WITH ZERO TO THE RIGHT

for D in [D1,D2,D3,D4,D5]:
    D['Ro_mid'] = np.hstack((D['Ro_mid'],10**1.5))
    D['RoH'] = np.hstack((D['RoH'],0.))

## CDISS
Ro = np.linspace(10**-3.5,10**1.5,10000)
cd = lambda d: 1/(1+(Ro/d))

## PLOT
fig,(ax1,ax2) = plt.subplots(2,1,figsize=(8,4.5),sharex=True)

s = 1e6

ax1.plot(D1['Ro_mid'],D1['RoH']/s,'C0',label=r'Low resolution, $\Delta x = $30km',lw=3)
ax1.plot(D2['Ro_mid'],D2['RoH']/16/s,'C2',label=r'High resolution, $\Delta x = $7.5km',lw=3)
ax1.plot(D3['Ro_mid'],D3['RoH']/s,'C1',label=r'LR + weak backscatter',ls='--')
ax1.plot(D4['Ro_mid'],D4['RoH']/s,'C3',label=r'LR + moderate backscatter',ls='--')
ax1.plot(D5['Ro_mid'],D5['RoH']/s,'C5',label=r'LR + strong backscatter',ls='--')

ytikmax = 0.25

ax1.axvline(np.log10(D1['Ro_mean']),0,ytikmax,c='C0',ls='-',lw=2)
ax1.axvline(np.log10(D2['Ro_mean']),0,ytikmax,c='C2',ls='-',lw=2)
ax1.axvline(np.log10(D3['Ro_mean']),0,ytikmax,c='C1',ls='--')
ax1.axvline(np.log10(D4['Ro_mean']),0,ytikmax,c='C3',ls='--')
ax1.axvline(np.log10(D5['Ro_mean']),0,ytikmax,c='C5',ls='--')

ax1.text(np.log10(D1['Ro_mean'])-0.01,1.3,'mean($R_o$)',rotation=90,ha='right',color='k')

ax2.plot(np.log10(Ro),cd(1.),'C1--',label=r'$c_{diss}$ for $R_{diss} = 1$, weak backscatter')
ax2.plot(np.log10(Ro),cd(6.),'C3--',label=r'$c_{diss}$ for $R_{diss} = 6$, moderate backscatter')
ax2.plot(np.log10(Ro),np.ones_like(Ro),'C5--',label=r'$c_{diss}$ for $R_{diss} = \infty$, strong backscatter')

ax2.plot(np.ones(2)*np.log10(1),[-1,0.5],'C1')
ax2.plot(np.ones(2)*np.log10(6),[-1,0.5],'C3')
ax2.text(np.log10(1)-0.03,0.3,r'$R_{diss} = 1$',rotation=90,ha='right',color='k')
ax2.text(np.log10(6)-0.03,0.3,r'$R_{diss} = 6$',rotation=90,ha='right',color='k')
ax2.text(-2.9,0.6,r'$c_{diss} = (1 + \frac{R_o}{R_{diss}})^{-1}$', fontsize=15)

ax2.legend(loc=3)
ax1.legend(loc=1)

ax1.set_xlim(-3,1)
ax1.set_ylim(-0.1,3.5)
ax2.set_ylim(-0.05,1.05)

ax1.set_title('a',loc='left',fontweight='bold')
ax1.set_title('Rossby number histogram')

ax2.set_title('b',loc='left',fontweight='bold')
ax2.set_title('Rossby number dissipation scaling')

ax2.set_xlabel('log$_{10}(R_o)$')
ax1.set_ylabel(r'$N$ $[10^6]$')

plt.tight_layout()
plt.savefig(path+'plots/Ro_hist.pdf')
plt.close(fig)
