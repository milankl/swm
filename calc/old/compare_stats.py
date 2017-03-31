## COMPARE STATISTICS

import numpy as np
import matplotlib.pyplot as plt

# load stats
stats1 = np.load('/home/mkloewer/python/swm/data/run0083/stats.npy').all()
stats2 = np.load('/home/mkloewer/python/swm/data/run0084/stats.npy').all()

time = np.array(stats1['time']) / 3600. / 24. # in days

key = 'PE'
quant1 = np.array(stats1[key]) #/ stats1[key][0]
quant2 = np.array(stats2[key]) #/ stats2[key][0] 

fig,(ax1,ax2) = plt.subplots(2,1,sharex=True)

ax1.plot(time,quant1,'b')
ax1.plot(time,quant2,'g')

ax1.set_ylabel('relative change')
ax1.set_title('Kinetic Energy')

key = 'KE'
quant1 = np.array(stats1[key]) #/ stats1[key][0]
quant2 = np.array(stats2[key]) #/ stats2[key][0] 

ax2.plot(time,quant1,'b',label=r'$B\nabla^4\mathbf{u}$')
ax2.plot(time,quant2,'g',label=r'$B[h^{-1}\nabla \cdot h\nabla]^2\mathbf{u}$')

ax2.set_xlabel('time [days]')
ax2.set_ylabel('relative change')
ax2.set_title('Potential Enstrophy')

ax2.legend(loc=1)

plt.tight_layout()
plt.show()

