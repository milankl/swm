## RIGHT HAND SIDE OF THE EQUATIONS
def rhs(u,v,h):
    """ Set of equations:
    
    u_t = qhv - P_x + Fx + Mx(u,v)
    v_t = -qhu - P_y + My(u,v)
    h_t = -(uh)_x - (vh)_y
    
    with P = .5*(u**2 + v**2) + gh, the bernoulli potential
    and q = (v_x - u_y + f)/h the potential vorticity
    
    using the enstrophy and energy conserving scheme (Arakawa and Lamb, 1981) and
    a friction term from Shchepetkin and O'Brien (1996) including their 
    increased accuracy of the derivative at the boundary.
    """
    
    h_u = ITu.dot(h)    # h on u-grid
    h_v = ITv.dot(h)    # h on v-grid    
    h_q = ITq.dot(h)    # h on q-grid

    U = u*h_u    # volume fluxes: U on u-grid
    V = v*h_v    # and V on v-grid
    
    q = (f_q + Gvx.dot(v) - Guy.dot(u)) / h_q   # potential vorticity q
    P = .5*(IuT.dot(u**2) + IvT.dot(v**2)) + g*h  # Bernoulli potential P

    # simple bi-harmonic mixing
    #diff_u = param['B']*LLu.dot(u)
    #diff_v = param['B']*LLv.dot(v)

    diff_u, diff_v = mixing(u,v,h,h_q,h_u,h_v)

    # Sadourny, 1975 enstrophy conserving scheme.
    #adv_u = Iqu.dot(q)*Ivu.dot(V)
    #adv_v = -Iqv.dot(q)*Iuv.dot(U)

    adv_u, adv_v = ALadvection(q,U,V)
    
    # add the terms to the right hand side
    rhs_u = adv_u - GTx.dot(P) + Fx + diff_u
    rhs_v = adv_v - GTy.dot(P) + diff_v
    rhs_h = -Gux.dot(U) - Gvy.dot(V)
    
    return rhs_u, rhs_v, rhs_h

def ALadvection(q,U,V):
    """ Arakawa and Lamb,1981 advection terms. See interpolation.py for further information. """
    
    AL1q = AL1.dot(q)
    AL2q = AL2.dot(q)
    
    adv_u = Seur.dot(ALeur.dot(q)*U) + Seul.dot(ALeul.dot(q)*U) +\
            Sau.dot(AL1q[indx_au]*V) + Sbu.dot(AL2q[indx_bu]*V) +\
            Scu.dot(AL2q[indx_cu]*V) + Sdu.dot(AL1q[indx_du]*V)

    adv_v = Spvu.dot(ALpvu.dot(q)*V) + Spvd.dot(ALpvd.dot(q)*V) -\
            Sav.dot(AL1q[indx_av]*U) - Sbv.dot(AL2q[indx_bv]*U) -\
            Scv.dot(AL2q[indx_cv]*U) - Sdv.dot(AL1q[indx_dv]*U)
    
    return adv_u, adv_v

def mixing(u,v,h,h_q,h_u,h_v):
    """ Biharmonic operator following the ideas of the divergence of a tensor approach by Shchepetkin and O'Brien, Mon. Wea. Rev. 1996. """
    # increased accuracy of the stencil at the boundary by using G2vx, G2uy
    
    P11 = h*(Gux.dot(u) - Gvy.dot(v))   # = -P22, component on the diagonal
    P12 = h_q*(G2vx.dot(v) + G2uy.dot(u))   # = P21, component on the off-diagonal
    diffu = (GTx.dot(P11) + Gqy.dot(P12)) / h_u
    diffv = (Gqx.dot(P12) - GTy.dot(P11)) / h_v
    
    # apply operator twice
    P11 = h*(Gux.dot(diffu) - Gvy.dot(diffv))   # = -P22, component on the diagonal
    P12 = h_q*(G2vx.dot(diffv) + G2uy.dot(diffu))   # = P21, component on the off-diagonal
    diff_u = param['B']*(GTx.dot(P11) + Gqy.dot(P12)) / h_u
    diff_v = param['B']*(Gqx.dot(P12) - GTy.dot(P11)) / h_v
    
    return diff_u, diff_v
