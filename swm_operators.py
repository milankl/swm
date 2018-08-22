## FUNCTIONs TO SET UP ALL OPERATORS

# TODO setting up of operators more readable!
# TODO set up some operators without 1/dx,1/dy to avoid division by large number
# TODO and then multiplication of it again

def set_grad_mat():
    """ Sets up the gradient matrices and makes them available as global variables.
    All matrices are set up as column sparse row.

    Notation:
        G           -   gradient
        T,u,v,q     -   on the T,u,v,q-grid
        x,y         -   in x- or y-direction.

    Hence, GTx is the x-derivative on the T-grid.

    See the shallow water model documentation for further information.
    """

    global GTx,GTy,Gux,Guy,Gvx,Gvy,Gqy,Gqx
    global G2vx,G2uy  # as above but with higher order stencils on the boundary

    # for readability unpack the param dictionary
    dx,dy = param['dx'],param['dy']
    nx,ny = param['nx'],param['ny']
    NT,Nu,Nv,Nq = param['NT'],param['Nu'],param['Nv'],param['Nq']

    # tangential velocity boundary condition
    lbc = param['lbc']

    # index used to delete the rows that correspond to a derivate across the east-west boundaries, i.e. to remove the periodicity in x
    indx1 = list(range(NT)); del indx1[(nx-1)::nx]
    indx2 = list(range(Nv+(ny-1))); del indx2[nx::(nx+1)]

    # d/dx from T-points to u-points
    GTx = (sparse.diags(np.ones(NT-1),1) - sparse.eye(NT))[indx1,:] / dx

    # d/dy from T-points to v-points
    GTy = (sparse.diags(np.ones(Nv),nx) - sparse.eye(NT))[:-nx,:] / dy

    # d/dx from u-points onto T-points
    Gux = -GTx.T.tocsr()

    # d/dy from v-points onto T-points
    Gvy = -GTy.T.tocsr()

    # d/dy from u-points to q-points
    # du/dy is zero at the western and eastern boundary
    d1 = np.ones(Nq) # set up the diagonal (including BC)
    d1[::(nx+1)] = 0. # western boundary
    d1[nx::(nx+1)] = 0. # eastern boundary
    indx3 = (d1 != 0) # index to remove unnecessary columns
    d1[-nx:-1] = lbc    # the boundary condition (north and south)
    Guy1 = sparse.diags(d1,0).tocsr()[:,indx3][:,nx-1:]
    Guy2 = sparse.diags(d1[::-1],0).tocsr()[:,indx3][:,:-(nx-1)] # fliplr and flipud of Guy1
    Guy = (Guy2 - Guy1) / dy

    # d/dy from q-points to u-points
    d1[-nx:-1] = 1.
    Gqy1 = sparse.diags(d1,0).tocsr()[:,indx3][:,nx-1:]
    Gqy2 = sparse.diags(d1[::-1],0).tocsr()[:,indx3][:,:-(nx-1)] # fliplr and flipud of Guy1
    Gqy = ( Gqy1 - Gqy2).T.tocsr()/ dy

    # d/dx from v-points to q-points
    # dv/dx is zero at the northern and southern boundary
    sj = Nv+(ny-1) # shape of Gvx in j-direction
    d2 = np.ones(sj) # set up the diagonal
    d2[::(nx+1)] = lbc  # boundary condition (east and west)
    Gvx = (sparse.dia_matrix((d2,-(nx+1)),shape=((Nq,sj))).tocsr()[:,indx2] +\
    sparse.dia_matrix((-d2[::-1],-(nx+1)),shape=((Nq,sj))).tocsr()[:,-np.array(indx2)[::-1]-1])/dx

    # d/dx from q-points to v-points
    # make use of Gvx
    d2[::(nx+1)] = 1.
    Gqx = (sparse.dia_matrix((d2,-(nx+1)),shape=((Nq,sj))).tocsr()[:,indx2] +\
    sparse.dia_matrix((-d2[::-1],-(nx+1)),shape=((Nq,sj))).tocsr()[:,-np.array(indx2)[::-1]-1])/dx
    Gqx = -Gqx.T.tocsr()

    ## HIGHER ORDER GRADIENTS - 2nd ORDER also on the boundary
    # stencil taken from Shchepetkin and O'Brien (1996), Mon. Wea. Review, equation (18)
    lateral_stencil = [4, -1,1/5.]

    # G2vx the 2nd order (also on the boundaries) gradient from v- on q-points
    block1 = np.array(lateral_stencil) # corresponding to the derivative at westernmost q-point
    block2 = np.array([[-1,1],]*(nx-1)) # corresponding to the interior derivatives
    block3 = -block1[::-1] # corresponding to the derivative at easternmost q-point

    datablock = np.hstack((block1,block2.flatten(),block3))
    data = np.array([datablock,]*(ny-1)).flatten() # data array

    lb = len(block1)
    row_ind_1st_row = np.hstack(([0,]*lb,np.arange(1,nx).repeat(2),[nx,]*lb)) + (nx+1)
    row_ind = (np.array([row_ind_1st_row,]*(ny-1)).T + np.arange(ny-1)*(nx+1)).T.flatten()

    # column indices
    col_ind_1st_row = np.hstack((range(lb),np.arange(nx).repeat(2)[1:-1],nx-np.arange(lb,0,-1)))
    col_ind = (np.array([col_ind_1st_row,]*(ny-1)).T + np.arange(ny-1)*nx).T.flatten()

    G2vx = sparse.csr_matrix((data,(row_ind,col_ind)),shape=(Nq,Nv)) / dx

    # G2uy
    block1 = np.array(lateral_stencil*(nx-1))    # derivative at southern boundary
    block2 = np.array([-1,1]*(nx-1)*(ny-1))     # derivative at interior points
    block3 = -block1[::-1]      # derivative at northern boundary

    data = np.hstack((block1,block2,block3))

    lb = int(len(block1)/(nx-1))
    row_ind_1st_row = np.arange(1,nx).repeat(lb)
    row_ind_last_row = Nq - row_ind_1st_row[::-1] - 1
    row_ind_2nd_row = np.arange(1,nx).repeat(2)
    row_ind_int = (np.array([row_ind_2nd_row,]*(ny-1)).T + np.arange(1,ny)*(nx+1)).T.flatten()
    row_ind = np.hstack((row_ind_1st_row,row_ind_int,row_ind_last_row))

    col_ind_1st_row = np.arange((nx-1)*lb).reshape(lb,nx-1).flatten('F')
    col_ind_int = (np.array([np.arange((nx-1)*(ny-1)),]*2).T + np.array([0,nx-1])).flatten()
    col_ind_last_row = Nu - col_ind_1st_row[::-1] - 1
    col_ind = np.hstack((col_ind_1st_row,col_ind_int,col_ind_last_row))

    G2uy = sparse.csr_matrix((data,(row_ind,col_ind)),shape=(Nq,Nu)) / dy

    # set data type
    Gux = Gux.astype(param['dat_type'])
    Guy = Guy.astype(param['dat_type'])
    Gvx = Gvx.astype(param['dat_type'])
    Gvy = Gvy.astype(param['dat_type'])
    GTx = GTx.astype(param['dat_type'])
    GTy = GTy.astype(param['dat_type'])
    Gqx = Gqx.astype(param['dat_type'])
    Gqy = Gqy.astype(param['dat_type'])
    G2uy = G2uy.astype(param['dat_type'])
    G2vx = G2vx.astype(param['dat_type'])

    if lbc != 2:    # higher order stencil is only for no-slip
        G2uy = Guy
        G2vx = Gvx


def set_lapl_mat():
    """ Sets up the horizontal Laplacian (harmonic diffusion) and also the
    biharmonic diffusion operator LL.

    Let Dx = d/dx, Dy = d/dy then

        Lu = Dx**2 + Dy**2 applied on the u-grid

    similar for Lv,LT,Lq.

    LLu = Lu*Lu = Dx**4 + Dy**4 + 2*Dx**2*Dy**2

    similar for LLv, LLT. The lateral boundary conditions are taken into account
    as they are contained in the G-matrices."""

    global Lu, Lv, Lq, LT
    global LLu, LLv, LLT, LLq

    # harmonic operators
    Lu = GTx.dot(Gux) + Gqy.dot(G2uy)
    Lv = Gqx.dot(G2vx) + GTy.dot(Gvy)
    LT = Gux.dot(GTx) + Gvy.dot(GTy)
    Lq = G2vx.dot(Gqx) + G2uy.dot(Gqy)

    # biharmonic operators
    LLu = Lu.dot(Lu)
    LLv = Lv.dot(Lv)
    LLT = LT.dot(LT)
    LLq = Lq.dot(Lq)

## SET UP INTERPOLATION MATRICES ON ARAKAWA C-grid
def set_interp_mat():
    """Sets up all 2- or 4-point interpolation between the u-,v-,T- and q-grid
    and makes them available as global matrices.

    Notation:
        I                           - Interpolation
        1. subscript u,v,q,T        - from grid u,v,q,T
        2. subscript u,v,q,T        - to grid u,v,q,T

    Hence, the Ivu interpolation matrix interpolates a variable from the v- onto
    the u-grid.

    """

    global Ivu, Iuv, IqT, IuT, IvT, ITu, ITv, Iqu, Iqv, Iuq, Ivq, ITq

    # for readability unpack the param dictionary
    dx,dy = param['dx'],param['dy']
    nx,ny = param['nx'],param['ny']
    NT,Nu,Nv,Nq = param['NT'],param['Nu'],param['Nv'],param['Nq']
    lbc = param['lbc']

    # index used to delete the rows that correspond to an interpolation across the east-west boundaries, i.e. to remove the periodicity in x
    indx1 = list(range(Nv+nx)); del indx1[(nx-1)::nx]
    indx2 = list(range(NT+ny-1)); del indx2[nx::(nx+1)]

    # interpolate v-points to u-points, 4-point average
    # including the information of the kinematic boundary condition
    d = np.ones(Nv) / 4. # diagonal
    Ivu = (sparse.dia_matrix((d,0),shape=(Nv+nx,Nv)) + \
    sparse.dia_matrix((d,1),shape=(Nv+nx,Nv)) +\
    sparse.dia_matrix((d,-nx+1),shape=(Nv+nx,Nv)) +\
    sparse.dia_matrix((d,-nx),shape=(Nv+nx,Nv)))[indx1,:]

    # interpolate u-points to v-points, 4-point average
    # including the information of the kinematic boundary condition
    Iuv = Ivu.T

    # interpolate q-points to T-points, 4-point average
    d = np.ones(Nq+ny-1) / 4. #diagonal
    IqT = (sparse.dia_matrix((d,0),shape=((NT+ny-1,Nq))) +\
    sparse.dia_matrix((d,1),shape=((NT+ny-1,Nq))) +\
    sparse.dia_matrix((d,nx+1),shape=((NT+ny-1,Nq))) +\
    sparse.dia_matrix((d,nx+2),shape=((NT+ny-1,Nq))))[indx2,:]

    # interpolate u-points to T-points, 2-point average
    IuT = abs(Gux*dx/2.)

    # interpolate v-points to T-points, 2-point average
    IvT = abs(Gvy*dy/2.)

    # interpolate T-points to u-points, 2-point average
    ITu = abs(GTx*dx/2.)

    # interpolate T-points to v-points, 2-point average
    ITv = abs(GTy*dy/2.)

    # interpolate q-points to u-points, 2-point average
    d = np.ones(Nq)/2.
    indx3 = list(range(Nq-nx-1)); del indx3[::(nx+1)]; del indx3[nx-1::nx]
    Iqu = (sparse.dia_matrix((d,0),shape=((Nq-nx-1,Nq))) +\
    sparse.dia_matrix((d,nx+1),shape=((Nq-nx-1,Nq))))[indx3,:]

    # interpolate u-points to q-points, 2-point average
    # include lateral boundary condition information in Iqu.T
    Iuq = Iqu.T.tocsr().copy()
    Iuq.data[:nx-1] = 1-lbc/2.
    Iuq.data[-nx+1:] = 1-lbc/2.

    # interpolate q-points to v-points, 2-point average
    # same diagonal d as for Iqu, reuse
    indx4 = list(range(Nv+ny-2)); del indx4[nx::(nx+1)]
    Iqv = (sparse.dia_matrix((d,nx+1),shape=((Nv+ny-2,Nq))) +\
    sparse.dia_matrix((d,nx+2),shape=((Nv+ny-2,Nq))))[indx4,:]

    # interpolate v-points to q-points, 2-point average
    # include lateral boundary condition information in Iqv.T
    Ivq = Iqv.T.tocsr().copy()
    Ivq.data[::2*nx] = 1-lbc/2.
    Ivq.data[::-2*nx] = 1-lbc/2.

    # interpolate T-points to q-points, copy T points to ghost points (no h gradients across boundaries)
    # data vector with entries increased by *2,*4 for ghost-point copy
    # equivalently: ITq.sum(axis=1) must be 1 in each row.
    d = np.ones(4*NT) / 4.
    d[1:2*nx-1] = .5
    d[-2*nx+1:-1] = .5
    d[2*nx:-2*nx-1:4*nx] = .5
    d[2*nx+1:-2*nx-1:4*nx] = .5
    d[-2*nx-1:2*nx:-4*nx] = .5
    d[-2*nx-2:2*nx:-4*nx] = .5
    d[[0,2*nx-1,-(2*nx),-1]] = 1

    ITq = sparse.csr_matrix(IqT.T)
    ITq.data = d

    # set data type, single/double precision
    IuT = IuT.astype(param['dat_type'])
    Iuv = Iuv.astype(param['dat_type'])
    Iuq = Iuq.astype(param['dat_type'])
    IvT = IvT.astype(param['dat_type'])
    Ivu = Ivu.astype(param['dat_type'])
    Ivq = Ivq.astype(param['dat_type'])
    ITu = ITu.astype(param['dat_type'])
    ITv = ITv.astype(param['dat_type'])
    ITq = ITq.astype(param['dat_type'])
    IqT = IqT.astype(param['dat_type'])
    Iqu = Iqu.astype(param['dat_type'])
    Iqv = Iqv.astype(param['dat_type'])

## RESHAPE FUNCTIONS
""" reshape functions are used to get from vector representation to matrix
representation mostly for plotting or output purposes. Note on the grid numbering:

    Any variable is numbered row-first.

    A_ij = (a b     =>      a = (a,b,c,d)
            c d)

    where physically i represents the y-axis and j the x-axis."""

def h2mat(eta):
    return eta.reshape((param['ny'],param['nx']))

def u2mat(u):
    return u.reshape((param['ny'],param['nx']-1))

def v2mat(v):
    return v.reshape((param['ny']-1,param['nx']))

def q2mat(q):
    return q.reshape((param['ny']+1,param['nx']+1))

## ARAKAWA and LAMB 1981 matrices
def set_arakawa_mat():
    """ Set up the linear combinations of potential vorticity as in
        Arakawa and Lamb 1981, eq. (3.34) as matrices. For the U-comp the arrangment is as follows

        c   a
     el   u   er
        d   b

        and for the V-component

          pu
        b   a
          v
        d   c
          pd

        in physical x,y space.

        """

    #TODO ALeur, ALeul, ALpvu, ALpvd can also be summarized in essentially two
    #TODO stencils, AL3, AL4, might lead to some minor speed up (see documentation)

    #TODO think about whether the S-matrices can also be written as index
    #TODO may involve the need of another ghost point though

    global AL1, AL2 # general matrices for the only two different stencils
    global indx_au, indx_bu, indx_cu, indx_du   # indices to pick the right rows of AL1 or 2
    global indx_av, indx_bv, indx_cv, indx_dv

    global ALeur, ALeul #U-components lin comb of q
    global Seul, Seur, Sau, Sbu, Scu, Sdu   #U-components shift operators

    global ALpvu, ALpvd #V-components lin comb of q
    global Spvu, Spvd, Sav, Sbv, Scv, Sdv   #V-components shift operators

    # for readability unpack the param dictionary
    nx,ny = param['nx'],param['ny']
    NT,Nu,Nv,Nq = param['NT'],param['Nu'],param['Nv'],param['Nq']

    ## FIRST AL1 and SECOND AL2 linear combination
    # AL1 corresponds to the 2 1 stencil, AL2 to the 1 2 stencil
    #                        1 2                     2 1
    d = np.ones(Nq+ny-1) / 24.  # data vector
    indx1 = list(range(Nq-nx-1)); del indx1[nx::(nx+1)]
    AL1 = (sparse.dia_matrix((2*d,0),shape=(Nq+nx,Nq)) +\
        sparse.dia_matrix((d,1),shape=(Nq+nx,Nq)) +\
        sparse.dia_matrix((d,nx+1),shape=(Nq+nx,Nq)) +\
        sparse.dia_matrix((2*d,nx+2),shape=(Nq+nx,Nq)))[indx1,:]

    # indices to pick from AL1/2 the correct rows for the a,b,c,d linear combination
    # for both U and V components
    indx_au = slice(-nx)
    indx_du = slice(nx,None)
    indx_av = list(range(nx**2)); del indx_av[nx-1::nx]
    indx_dv = list(range(1,nx**2+1)); del indx_dv[nx-1::nx]
    indx_av = np.array(indx_av)     # convert to numpy array for speed up
    indx_dv = np.array(indx_dv)

    AL2 = (sparse.dia_matrix((d,0),shape=(Nq+nx,Nq)) +\
        sparse.dia_matrix((2*d,1),shape=(Nq+nx,Nq)) +\
        sparse.dia_matrix((2*d,nx+1),shape=(Nq+nx,Nq)) +\
        sparse.dia_matrix((d,nx+2),shape=(Nq+nx,Nq)))[indx1,:]

    indx_bu = indx_du
    indx_cu = indx_au
    indx_bv = indx_dv
    indx_cv = indx_av

    ## U-component of advection
    # ALeur is the epsilon linear combination to the right of the associated u-point
    indx2 = list(range(NT+ny-1)); del indx2[nx::(nx+1)]; del indx2[nx-1::nx]
    ALeur = (sparse.dia_matrix((d,0),shape=(NT+ny-1,Nq)) +\
        sparse.dia_matrix((d,1),shape=(NT+ny-1,Nq)) +\
        sparse.dia_matrix((-d,nx+1),shape=(NT+ny-1,Nq)) +\
        sparse.dia_matrix((-d,nx+2),shape=(NT+ny-1,Nq)))[indx2,:]

    # ALeul is the epsilon linear combination to the left of the associated u-point
    indx3 = list(range(NT+ny-1)); del indx3[nx::(nx+1)]; del indx3[::nx]
    ALeul = (sparse.dia_matrix((-d,0),shape=(NT+ny-1,Nq)) +\
        sparse.dia_matrix((-d,1),shape=(NT+ny-1,Nq)) +\
        sparse.dia_matrix((d,nx+1),shape=(NT+ny-1,Nq)) +\
        sparse.dia_matrix((d,nx+2),shape=(NT+ny-1,Nq)))[indx3,:]

    # Seur, Seul are shift-matrices so that the correct epsilon term is taken for the associated u-point
    ones = np.ones(Nu)
    ones[::(nx-1)] = 0
    Seur = sparse.dia_matrix((ones,1),shape=(Nu,Nu)).tocsr()
    Seul = Seur.T.tocsr()

    # Shift matrices for the a,b,c,d linear combinations
    ones = np.ones(Nu+ny)
    indx4 = list(range(Nv+nx)); del indx4[(nx-1)::nx]

    Sau = sparse.dia_matrix((ones,1),shape=(Nu+ny,Nv)).tocsr()[indx4,:]
    Scu = sparse.dia_matrix((ones,0),shape=(Nu+ny,Nv)).tocsr()[indx4,:]
    Sbu = sparse.dia_matrix((ones,-(nx-1)),shape=(Nu+ny,Nv)).tocsr()[indx4,:]
    Sdu = sparse.dia_matrix((ones,-nx),shape=(Nu+ny,Nv)).tocsr()[indx4,:]


    ## V-component of advection
    # ALpvu is the p linear combination, up
    indx5 = list(range(Nq-nx-1)); del indx5[nx::(nx+1)]; del indx5[-nx:]
    ALpvu = (sparse.dia_matrix((-d,0),shape=(Nq,Nq)) +\
        sparse.dia_matrix((d,1),shape=(Nq,Nq)) +\
        sparse.dia_matrix((-d,nx+1),shape=(Nq,Nq)) +\
        sparse.dia_matrix((d,nx+2),shape=(Nq,Nq)))[indx5,:]

    # ALpvd is the p linear combination, down
    indx6 = list(range(Nq-nx-1)); del indx6[nx::(nx+1)]; del indx6[:nx]
    ALpvd = (sparse.dia_matrix((d,0),shape=(Nq,Nq)) +\
        sparse.dia_matrix((-d,1),shape=(Nq,Nq)) +\
        sparse.dia_matrix((d,nx+1),shape=(Nq,Nq)) +\
        sparse.dia_matrix((-d,nx+2),shape=(Nq,Nq)))[indx6,:]

    # associated shift matrix
    ones = np.ones(Nv)
    Spvu = sparse.dia_matrix((ones,nx),shape=(Nv,Nv)).tocsr()
    Spvd = Spvu.T.tocsr()

    # Shift matrices for a,b,c,d linear combinations
    ones = np.ones(Nv+nx)
    indx7 = list(range(Nv+nx)); del indx7[(nx-1)::nx]

    Sav = sparse.dia_matrix((ones,-nx),shape=(Nv+nx,Nv)).tocsr()[indx7,:].T.tocsr()
    Sbv = sparse.dia_matrix((ones,-(nx-1)),shape=(Nv+nx,Nv)).tocsr()[indx7,:].T.tocsr()
    Scv = sparse.dia_matrix((ones,0),shape=(Nv+nx,Nv)).tocsr()[indx7,:].T.tocsr()
    Sdv = sparse.dia_matrix((ones,1),shape=(Nv+nx,Nv)).tocsr()[indx7,:].T.tocsr()

    # set data type single/double precision
    AL1 = AL1.astype(param['dat_type'])
    AL2 = AL2.astype(param['dat_type'])

    ALeur = ALeur.astype(param['dat_type'])
    ALeul = ALeul.astype(param['dat_type'])
    Seul = Seul.astype(param['dat_type'])
    Seur = Seur.astype(param['dat_type'])
    Sau = Sau.astype(param['dat_type'])
    Sbu = Sbu.astype(param['dat_type'])
    Scu = Scu.astype(param['dat_type'])
    Sdu = Sdu.astype(param['dat_type'])

    ALpvu = ALpvu.astype(param['dat_type'])
    ALpvd = ALpvd.astype(param['dat_type'])
    Spvu = Spvu.astype(param['dat_type'])
    Spvd = Spvd.astype(param['dat_type'])
    Sav = Sav.astype(param['dat_type'])
    Sbv = Sbv.astype(param['dat_type'])
    Scv = Scv.astype(param['dat_type'])
    Sdv = Sdv.astype(param['dat_type'])
