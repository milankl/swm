## HISTOGRAM PLOTTING FOR REYNOLDS AND ROSSBY NUMBERS
from __future__ import print_function
path = '/home/mkloewer/python/swm/'
import os; os.chdir(path) # change working directory
import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt
import time as tictoc
from netCDF4 import Dataset
import glob
import matplotlib

plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'

runfolder = [0,6,10,14,15]  # with bottom friction folders
print('Creating ReRo histogram plot from run '+str(runfolder)) 

# read data
runpath = path+'data/run%04i' % runfolder[0]
D1 = np.load(runpath+'/analysis/Re_hist.npy').all()   
D1.update(np.load(runpath+'/analysis/Ro_hist.npy').all())

runpath = path+'data/run%04i' % runfolder[1]
D2 = np.load(runpath+'/analysis/Re_hist.npy').all()
D2.update(np.load(runpath+'/analysis/Ro_hist.npy').all())

runpath = path+'data/run%04i' % runfolder[2]
D3 = np.load(runpath+'/analysis/Re_hist.npy').all()
D3.update(np.load(runpath+'/analysis/Ro_hist.npy').all())

runpath = path+'data/run%04i' % runfolder[3]
D4 = np.load(runpath+'/analysis/Re_hist.npy').all()
D4.update(np.load(runpath+'/analysis/Ro_hist.npy').all())

runpath = path+'data/run%04i' % runfolder[4]
D5 = np.load(runpath+'/analysis/Re_hist.npy').all()
D5.update(np.load(runpath+'/analysis/Ro_hist.npy').all())

# OPTIONS
runfolder = [3,10,13,12,14] #without bottom friction folders

## read data
runpath = path+'data/newold/run%04i' % runfolder[0]
D11 = np.load(runpath+'/analysis/Re_hist.npy').all()   
D11.update(np.load(runpath+'/analysis/Ro_hist.npy').all())

runpath = path+'data/newold/run%04i' % runfolder[1]
D12 = np.load(runpath+'/analysis/Re_hist.npy').all()
D12.update(np.load(runpath+'/analysis/Ro_hist.npy').all())

runpath = path+'stoch/data/run%04i' % runfolder[2]
D13 = np.load(runpath+'/analysis/Re_hist.npy').all()
D13.update(np.load(runpath+'/analysis/Ro_hist.npy').all())

runpath = path+'stoch/data/run%04i' % runfolder[3]
D14 = np.load(runpath+'/analysis/Re_hist.npy').all()
D14.update(np.load(runpath+'/analysis/Ro_hist.npy').all())

runpath = path+'stoch/data/run%04i' % runfolder[4]
D15 = np.load(runpath+'/analysis/Re_hist.npy').all()
D15.update(np.load(runpath+'/analysis/Ro_hist.npy').all())

## PLOT
fig,axs = plt.subplots(2,2,figsize=(12,8))

axs[0,0].plot(D1['Ro_mid'],D1['RoH'],'C0',lw=3)
axs[0,0].plot(D2['Ro_mid'],D2['RoH']/16,'C2',lw=3)
axs[0,0].plot(D3['Ro_mid'],D3['RoH'],'C3',ls='--')
axs[0,0].plot(D4['Ro_mid'],D4['RoH'],'C1',ls='--')
axs[0,0].plot(D5['Ro_mid'],D5['RoH'],'C5',ls='--')

l01, = axs[0,1].plot(D1['Re_mid'],D1['ReH'],'C0',label=r'Low resolution, $\Delta x = $30km',lw=3)
l02, = axs[0,1].plot(D2['Re_mid'],D2['ReH']/16,'C2',label=r'High resolution, $\Delta x = $7.5km',lw=3)
l03, = axs[0,1].plot(D3['Re_mid'],D3['ReH'],'C3',label=r'LR + weak backscatter',ls='--')
l04, = axs[0,1].plot(D4['Re_mid'],D4['ReH'],'C1',label=r'LR + moderate backscatter',ls='--')
l05, = axs[0,1].plot(D5['Re_mid'],D5['ReH'],'C5',label=r'LR + strong backscatter',ls='--')

axs[1,0].plot(D11['Ro_mid'],D11['RoH'],'C0',lw=3)
axs[1,0].plot(D12['Ro_mid'],D12['RoH']/16,'C2',lw=3)
axs[1,0].plot(D13['Ro_mid'],D13['RoH'],'C3',ls='--')
axs[1,0].plot(D14['Ro_mid'],D14['RoH'],'C1',ls='--')
axs[1,0].plot(D15['Ro_mid'],D15['RoH'],'C5',ls='--')

axs[1,1].plot(D11['Re_mid'],D11['ReH'],'C0',lw=3)
axs[1,1].plot(D12['Re_mid'],D12['ReH']/16,'C2',lw=3)
axs[1,1].plot(D13['Re_mid'],D13['ReH'],'C3',ls='--')
axs[1,1].plot(D14['Re_mid'],D14['ReH'],'C1',ls='--')
axs[1,1].plot(D15['Re_mid'],D15['ReH'],'C5',ls='--')

ytikmax = 0.08
ytikmin = 1-ytikmax

axs[0,0].axvline(np.log10(D1['Ro_mean']),0,ytikmax,c='C0',ls='-',lw=2,label=r'$\langle \overline{R_o} \rangle$ = %.3f' % D1['Ro_mean'])
axs[0,0].axvline(np.log10(D2['Ro_mean']),0,ytikmax,c='C2',ls='-',lw=2,label=r'$\langle \overline{R_o} \rangle$ = %.3f' % D2['Ro_mean'])
axs[0,0].axvline(np.log10(D3['Ro_mean']),0,ytikmax,c='C3',ls='--',label=r'$\langle \overline{R_o} \rangle$ = %.3f' % D3['Ro_mean'])
axs[0,0].axvline(np.log10(D4['Ro_mean']),0,ytikmax,c='C1',ls='--',label=r'$\langle \overline{R_o} \rangle$ = %.3f' % D4['Ro_mean'])
axs[0,0].axvline(np.log10(D5['Ro_mean']),0,ytikmax,c='C5',ls='--',label=r'$\langle \overline{R_o} \rangle$ = %.3f' % D5['Ro_mean'])

l1 = axs[0,1].axvline(np.log10(D1['Re_mean']),0,ytikmax,c='C0',ls='-',lw=2,label=r'$\langle \overline{R_e} \rangle$ = %i' % D1['Re_mean'])
l2 = axs[0,1].axvline(np.log10(D2['Re_mean']),0,ytikmax,c='C2',ls='-',lw=2,label=r'$\langle \overline{R_e} \rangle$ = %i' % D2['Re_mean'])
l3 = axs[0,1].axvline(np.log10(D3['Re_mean']),0,ytikmax,c='C3',ls='--',label=r'$\langle \overline{R_e} \rangle$ = %i' % D3['Re_mean'])
l4 = axs[0,1].axvline(np.log10(D4['Re_mean']),0,ytikmax,c='C1',ls='--',label=r'$\langle \overline{R_e} \rangle$ = %i' % D4['Re_mean'])
l5 = axs[0,1].axvline(np.log10(D5['Re_mean']),0,ytikmax,c='C5',ls='--',label=r'$\langle \overline{R_e} \rangle$ = %i' % D5['Re_mean'])

axs[1,0].axvline(np.log10(D11['Ro_mean']),0,ytikmax,c='C0',ls='-',lw=2,label=r'$\langle \overline{R_o} \rangle$ = %.3f' % D11['Ro_mean'])
axs[1,0].axvline(np.log10(D12['Ro_mean']),0,ytikmax,c='C2',ls='-',lw=2,label=r'$\langle \overline{R_o} \rangle$ = %.3f' % D12['Ro_mean'])
axs[1,0].axvline(np.log10(D13['Ro_mean']),0,ytikmax,c='C3',ls='--',label=r'$\langle \overline{R_o} \rangle$ = %.3f' % D13['Ro_mean'])
axs[1,0].axvline(np.log10(D14['Ro_mean']),0,ytikmax,c='C1',ls='--',label=r'$\langle \overline{R_o} \rangle$ = %.3f' % D14['Ro_mean'])
axs[1,0].axvline(np.log10(D15['Ro_mean']),0,ytikmax,c='C5',ls='--',label=r'$\langle \overline{R_o} \rangle$ = %.3f' % D15['Ro_mean'])

axs[1,1].axvline(np.log10(D11['Re_mean']),0,ytikmax,c='C0',ls='-',lw=2,label=r'$\langle \overline{R_e} \rangle$ = %i' % D11['Re_mean'])
axs[1,1].axvline(np.log10(D12['Re_mean']),0,ytikmax,c='C2',ls='-',lw=2,label=r'$\langle \overline{R_e} \rangle$ = %i' % D12['Re_mean'])
axs[1,1].axvline(np.log10(D13['Re_mean']),0,ytikmax,c='C3',ls='--',label=r'$\langle \overline{R_e} \rangle$ = %i' % D13['Re_mean'])
axs[1,1].axvline(np.log10(D14['Re_mean']),0,ytikmax,c='C1',ls='--',label=r'$\langle \overline{R_e} \rangle$ = %i' % D14['Re_mean'])
axs[1,1].axvline(np.log10(D15['Re_mean']),0,ytikmax,c='C5',ls='--',label=r'$\langle \overline{R_e} \rangle$ = %i' % D15['Re_mean'])

axs[0,0].axvline(np.log10(D1['Ro_mean']),ytikmin,1,c='C0',ls='-',lw=2)
axs[0,0].axvline(np.log10(D2['Ro_mean']),ytikmin,1,c='C2',ls='-',lw=2)
axs[0,0].axvline(np.log10(D3['Ro_mean']),ytikmin,1,c='C3',ls='--')
axs[0,0].axvline(np.log10(D4['Ro_mean']),ytikmin,1,c='C1',ls='--')
axs[0,0].axvline(np.log10(D5['Ro_mean']),ytikmin,1,c='C5',ls='--')

axs[0,1].axvline(np.log10(D1['Re_mean']),ytikmin,1,c='C0',ls='-',lw=2)
axs[0,1].axvline(np.log10(D2['Re_mean']),ytikmin,1,c='C2',ls='-',lw=2)
axs[0,1].axvline(np.log10(D3['Re_mean']),ytikmin,1,c='C3',ls='--')
axs[0,1].axvline(np.log10(D4['Re_mean']),ytikmin,1,c='C1',ls='--')
axs[0,1].axvline(np.log10(D5['Re_mean']),ytikmin,1,c='C5',ls='--')

actual_legend = plt.legend(handles=[l01,l02,l03,l04,l05],loc=2)
axs[1,1].add_artist(actual_legend)
axs[0,1].legend(handles=[l1,l2,l3,l4,l5],loc=1,fontsize=8,frameon=False)
axs[0,0].legend(loc=2,fontsize=8,frameon=False)
axs[1,0].legend(loc=2,fontsize=8,frameon=False)
axs[1,1].legend(loc=1,fontsize=8,frameon=False)

axs[1,0].axvline(np.log10(D11['Ro_mean']),ytikmin,1,c='C0',ls='-',lw=2)
axs[1,0].axvline(np.log10(D12['Ro_mean']),ytikmin,1,c='C2',ls='-',lw=2)
axs[1,0].axvline(np.log10(D13['Ro_mean']),ytikmin,1,c='C3',ls='--')
axs[1,0].axvline(np.log10(D14['Ro_mean']),ytikmin,1,c='C1',ls='--')
axs[1,0].axvline(np.log10(D15['Ro_mean']),ytikmin,1,c='C5',ls='--')

axs[1,1].axvline(np.log10(D11['Re_mean']),ytikmin,1,c='C0',ls='-',lw=2)
axs[1,1].axvline(np.log10(D12['Re_mean']),ytikmin,1,c='C2',ls='-',lw=2)
axs[1,1].axvline(np.log10(D13['Re_mean']),ytikmin,1,c='C3',ls='--')
axs[1,1].axvline(np.log10(D14['Re_mean']),ytikmin,1,c='C1',ls='--')
axs[1,1].axvline(np.log10(D15['Re_mean']),ytikmin,1,c='C5',ls='--')

axs[0,0].set_xlim(-3.5,0)
axs[0,0].set_ylim(1,3e6)

axs[1,0].set_xlim(-3.5,0)
axs[1,0].set_ylim(1,3e6)

axs[0,1].set_xlim(-2,4.5)
axs[0,1].set_ylim(1,3e6)

axs[1,1].set_xlim(-2,4.5)
axs[1,1].set_ylim(1,3e6)

axs[0,0].set_title('a',loc='left',fontweight='bold')
axs[0,1].set_title('b',loc='left',fontweight='bold')
axs[1,0].set_title('c',loc='left',fontweight='bold')
axs[1,1].set_title('d',loc='left',fontweight='bold')

axs[0,0].set_xticklabels([])
axs[0,1].set_xticklabels([])

axs[0,1].set_yticklabels([])
axs[1,1].set_yticklabels([])

axs[0,0].set_title('Rossby number histogram, with bottom friction')
axs[0,1].set_title('Reynolds number histogram, with bottom friction')
axs[1,0].set_title('Rossby number histogram, without bottom friction')
axs[1,1].set_title('Reynolds number histogram, without bottom friction')

axs[1,0].set_xlabel('log$_{10}(R_o)$')
axs[1,1].set_xlabel('log$_{10}(R_e)$')

axs[0,0].set_ylabel(r'$N\quad[10^6]$')
axs[1,0].set_ylabel(r'$N\quad[10^6]$')

axs[0,0].set_yticklabels(axs[0,0].get_yticks()/1e6)
axs[1,0].set_yticklabels(axs[1,0].get_yticks()/1e6)

plt.tight_layout()
plt.savefig(path+'compare/ReRo_hist.pdf')
plt.close(fig)