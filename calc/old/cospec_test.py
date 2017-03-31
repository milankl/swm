## test cospec

import numpy as np
import matplotlib.pyplot as plt


t = np.linspace(0,300,150)
x = np.linspace(0,200,112)

xx,tt = np.meshgrid(x,t)

dt = t[1] - t[0]
dx = x[1] - x[0]

k = -0.04
w = .015

xray = x[-1] + w/k*t

z = np.sin(2*np.pi*(k*xx - w*tt)) + 0*np.random.randn(len(t),len(x))

def cospec(A,dx,dt):
    """ assumed shape of A is time x space. """
    nt,nx = np.shape(A)
    
    f = np.fft.fftfreq(nt,dt)
    k = np.fft.fftfreq(nx,dx)
    
    idxf = int(np.ceil(nt/2.))    # index of maximum positive frequency
    idxk = int(np.ceil(nx/2.))    # index of maximum positive wavenumber

    # 2D FFT
    p = abs(np.fft.fft2(A))**2
    
    f = f[1:idxf]    # kill 0 and negative frequencies
    # reorder to have -kmax to +kmax
    k = np.hstack((k[idxk:],k[1:idxk]))
    # apply also for p
    p = np.hstack((p[1:idxf,idxk:],p[1:idxf,1:idxk]))[:,::-1]
    
    #TODO issue: p_max is shifted by -1 in positive k-direction

    return f,k,p

f,r,p = cospec(z,dx,dt)
i,j = np.unravel_index(np.argmax(p),p.shape)

fig,(ax1,ax2) = plt.subplots(1,2)

ax1.pcolormesh(x,t,z)
ax1.plot(xray,t,'w',lw=2)
ax1.set_xlim(x[0],x[-1])
ax1.set_ylim(t[0],t[-1])

ax2.pcolormesh(r,f,np.log10(p))
ax2.scatter(k,w,20,'w')
ax2.scatter(r[j],f[i],20,'r')
print(w,k)
print(f[i],r[j])

# ax2.set_xlim(-2,2)
# ax2.set_ylim(0,0.2)
plt.show()