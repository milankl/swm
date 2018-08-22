## Visualize Arakawa C-grid numbering
import numpy as np
import matplotlib.pyplot as plt

nx = 3
ny = 3

Lx = nx
Ly = ny

dx = Lx/nx
dy = Ly/ny

# grid vectors for T-points
x_T = np.arange(dx/2.,Lx,dx)
y_T = np.arange(dy/2.,Ly,dy)

# grid vectors for u-points
x_u = x_T[:-1] + dx/2.
y_u = y_T

#grid vectors for v-points
x_v = x_T
y_v = y_T[:-1] + dy/2.

# grid vectors for q-points
x_q = np.arange(0,Lx+dx/2.,dx)
y_q = np.arange(0,Ly+dy/2.,dy)

## meshgrid

xx_T,yy_T = np.meshgrid(x_T,y_T)
xx_u,yy_u = np.meshgrid(x_u,y_u)
xx_v,yy_v = np.meshgrid(x_v,y_v)
xx_q,yy_q = np.meshgrid(x_q,y_q)

xx_T = xx_T.flatten()
xx_u = xx_u.flatten()
xx_v = xx_v.flatten()
xx_q = xx_q.flatten()

yy_T = yy_T.flatten()
yy_u = yy_u.flatten()
yy_v = yy_v.flatten()
yy_q = yy_q.flatten()

fig,ax = plt.subplots(1,1)

ax.plot(xx_T,yy_T,'r.',ms=10,mec='k')
ax.plot(xx_u,yy_u,'g.',ms=10,mec='k')
ax.plot(xx_v,yy_v,'b.',ms=10,mec='k')
ax.plot(xx_q,yy_q,'k.',ms=10)

ax.plot(x_v,[0]*nx,'b.',ms=10,alpha=.6)
ax.plot(x_v,[Ly]*nx,'b.',ms=10,alpha=.6)
ax.plot([0]*ny,y_u,'g.',ms=10,alpha=.6)
ax.plot([Lx]*ny,y_u,'g.',ms=10,alpha=.6)


plt.grid()

for xx,yy,vv in zip([xx_T,xx_u,xx_v,xx_q],[yy_T,yy_u,yy_v,yy_q],['T','u','v','q']):
    for i in range(len(xx)):
        ax.text(xx[i]+0.03,yy[i]+0.03,r'$'+vv+'_{%i}$' %(i+1),fontsize=15)

ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')

ax.set_xticks(range(Lx+1))
ax.set_yticks(range(Ly+1))

ax.set_xlim(-0.1,nx+0.1)
ax.set_ylim(-0.1,ny+0.1)

ax.set_frame_on(False)

dxstring = [(r'$%i\Delta x$' % i) if i > 1 else r'$\Delta x$' for i in range(1,nx)]
dystring = [(r'$%i\Delta y$' % i) if i > 1 else r'$\Delta y$' for i in range(1,ny)]

ax.set_xticklabels([r'$0$']+dxstring+[r'$L_x$'])
ax.set_yticklabels([r'$0$']+dystring+[r'$L_y$'])
plt.tight_layout()
plt.show()