import numpy as np
import matplotlib.pyplot as plt

N = 10**7
x = np.random.randn(N) + 3.
y = np.random.randn(N) + 1.

##
H1,xe1,ye1 = np.histogram2d(x,y,bins=100,range=[[x.min(),x.max()],[y.min(),y.max()]],normed=True)
H2,xe2,ye2 = np.histogram2d(x,y,bins=100,range=[[x.min(),x.max()],[y.min()+3,y.max()]],normed=True)
H3,xe3,ye3 = np.histogram2d(x,y,bins=100,range=[[x.min(),x.max()],[y.min()-3,y.max()]],normed=True)

H4,ye4 = np.histogram(y,bins=100,range=[y.min(),y.max()],normed=True)
H5,ye5 = np.histogram(y,bins=100,range=[y.min()-3,y.max()],normed=True)
H6,ye6 = np.histogram(y,bins=100,range=[y.min()+3,y.max()],normed=True)

fig,(ax1,ax2,ax3) = plt.subplots(1,3,sharex=True,sharey=True)

c1 = ax1.pcolormesh(xe1,ye1,H1)
plt.colorbar(c1,ax=ax1)
ax1.set_xlim(1,5)
ax1.set_ylim(-1,3)

c2 = ax2.pcolormesh(xe2,ye2,H2)
plt.colorbar(c2,ax=ax2)

c3 = ax3.pcolormesh(xe3,ye3,H3)
plt.colorbar(c3,ax=ax3)

fig2,ax4 = plt.subplots(1,1)

ax4.plot(ye4[:-1],H4,drawstyle='steps')
ax4.plot(ye5[:-1],H5,drawstyle='steps')
ax4.plot(ye6[:-1],H6,drawstyle='steps')

plt.show()


