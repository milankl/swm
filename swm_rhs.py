## RIGHT HAND SIDE OF THE EQUATIONS
def rhs(u,v,eta,e):
    """ Set of equations:

    u_t = qhv - p_x + Fx + Mx(u,v) - bottom_friction
    v_t = -qhu - p_y + My(u,v)  - bottom_friction
    eta_t = -(uh)_x - (vh)_y

    with p = .5*(u**2 + v**2) + gh, the bernoulli potential
    and q = (v_x - u_y + f)/h the potential vorticity

    using the enstrophy and energy conserving scheme (Arakawa and Lamb, 1981) and
    a biharmonic lateral mixing term based on Shchepetkin and O'Brien (1996).
    """

    h = eta + param['H']

    h_u = ITu.dot(h)    # h on u-grid
    h_v = ITv.dot(h)    # h on v-grid
    h_q = ITq.dot(h)    # h on q-grid

    U = u*h_u    # volume fluxes: U on u-grid
    V = v*h_v    # and V on v-grid

    dudx = Gux.dot(u)   # precompute spatial derivatives
    dudy = Guy.dot(u)
    dvdx = Gvx.dot(v)
    dvdy = Gvy.dot(v)

    KE = IuT.dot(u**2) + IvT.dot(v**2)  # kinetic energy without .5-factor

    q = (f_q + dvdx - dudy) / h_q       # potential vorticity q
    p = .5*KE + param['g']*h            # Bernoulli potential p

    ## BOTTOM FRICTION: quadratic drag
    sqrtKE = np.sqrt(KE)
    bfric_u = param['c_D']*ITu.dot(sqrtKE)*u/h_u
    bfric_v = param['c_D']*ITv.dot(sqrtKE)*v/h_v

    ## ADVECTION
    # Sadourny, 1975 enstrophy conserving scheme.
    # adv_u = Iqu.dot(q)*Ivu.dot(V)
    # adv_v = -Iqv.dot(q)*Iuv.dot(U)

    # Arakawa and Lamb, 1981
    adv_u, adv_v = ALadvection(q,U,V)

    ## LATERAL MIXING OPERATOR
    # simple bi-harmonic mixing
    # bidiff_u = param['nu_B']*LLu.dot(u)
    # bidiff_v = param['nu_B']*LLv.dot(v)

    # symmetric stress tensor S = (S11, S12, S12, -S11), store only S11, S12
    S = (dudx - dvdy,G2vx.dot(v) + G2uy.dot(u))
    hS = (h*S[0],h_q*S[1])

    diff_u = (GTx.dot(hS[0]) + Gqy.dot(hS[1])) / h_u
    diff_v = (Gqx.dot(hS[1]) - GTy.dot(hS[0])) / h_v

    # biharmonic stress tensor R = (R11, R12, R12, -R11), store only R11, R12
    R = (Gux.dot(diff_u) - Gvy.dot(diff_v), G2vx.dot(diff_v) + G2uy.dot(diff_u))
    nuhR = (param['nu_B']*h*R[0],param['nu_B']*h_q*R[1])

    bidiff_u = (GTx.dot(nuhR[0]) + Gqy.dot(nuhR[1])) / h_u
    bidiff_v = (Gqx.dot(nuhR[1]) - GTy.dot(nuhR[0])) / h_v

    ## backscatter
    D2 = S[0]**2 + IqT.dot(S[1]**2)            # Deformation rate squared
    #c_diss = 1./(1+np.sqrt(D2)/f_T)**(param['n_diss'])   # D/f = Ro, the Rossby number for subgrid-EKE

    # NEW VERSION
    c_diss = 1./(1+(np.sqrt(D2)/f_T)**param['n_diss'])

    e_over_h = e/h
    nu_back = -param['c_back']*param['max_dxdy']*np.sqrt(2*e_over_h.clip(0,e_over_h.max()))
    #diff_e = param['nu_e']*(Gux.dot(h_u*GTx.dot(e)) + Gvy.dot(h_v*GTy.dot(e)))/h
    diff_e = param['nu_e']*LT.dot(e)
    E_back = nu_back*h*D2
    E_diss = c_diss*param['nu_B']*(hS[0]*R[0] + IqT.dot(hS[1]*R[1]))

    nu_back_hS0 = nu_back*hS[0]
    nu_back_hS1 = ITq.dot(nu_back)*hS[1]
    back_diff_u = (GTx.dot(nu_back_hS0) + Gqy.dot(nu_back_hS1)) / h_u
    back_diff_v = (Gqx.dot(nu_back_hS1) - GTy.dot(nu_back_hS0)) / h_v

    ## RIGHT-HAND SIDE: ADD TERMS
    rhs_u = adv_u - GTx.dot(p) + Fx/h_u - bidiff_u - bfric_u + back_diff_u
    rhs_v = adv_v - GTy.dot(p) - bidiff_v - bfric_v + back_diff_v
    rhs_eta = -(Gux.dot(U) + Gvy.dot(V))
    rhs_e = -E_diss + E_back + diff_e

    return rhs_u, rhs_v, rhs_eta, rhs_e

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
