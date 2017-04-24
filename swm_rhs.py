## RIGHT HAND SIDE OF THE EQUATIONS
def rhs(u,v,eta):
    """ Set of equations:
    
    u_t = qhv - p_x + Fx + Mx(u,v) - bottom_friction
    v_t = -qhu - p_y + My(u,v)  - bottom_friction
    eta_t = -(uh)_x - (vh)_y
    
    with p = .5*(u**2 + v**2) + gh, the bernoulli potential
    and q = (v_x - u_y + f)/h the potential vorticity
    
    using the enstrophy and energy conserving scheme (Arakawa and Lamb, 1981) and
    a biharmonic lateral mixing term based on Shchepetkin and O'Brien (1996).
    """
    
    #TODO param[nu_B] is large, applying the biharmonic creates tiny values (as dx^4 is large)
    #TODO think about a way to avoid possibly involved rounding errors especially with single precision
    #TODO might be efficiently only possible for dx=dy
    
    h = eta + param['H']    
    
    h_u = ITu*h    # h on u-grid
    h_v = ITv*h    # h on v-grid    
    h_q = ITq*h    # h on q-grid
    
    U = u*h_u    # volume fluxes: U on u-grid
    V = v*h_v    # and V on v-grid
    
    KE = IuT*(u**2) + IvT*(v**2)  # kinetic energy without .5-factor
    
    q = (f_q + Gvx*v - Guy*u) / h_q       # potential vorticity q
    p = .5*KE + param['g']*h            # Bernoulli potential p
    
    ## BOTTOM FRICTION: quadratic drag
    sqrtKE = np.sqrt(KE)
    bfric_u = param['c_D']*(ITu*sqrtKE)*u/h_u
    bfric_v = param['c_D']*(ITv*sqrtKE)*v/h_v
    
    ## ADVECTION
    # Sadourny, 1975 enstrophy conserving scheme.
    #adv_u = (Iqu*q)*(Ivu*V)
    #adv_v = -(Iqv*q)*(Iuv*U)
    
    # Arakawa and Lamb, 1981
    adv_u, adv_v = ALadvection(q,U,V)
    
    ## LATERAL MIXING OPERATOR
    # crude but simple bi-harmonic mixing
    #bidiff_u = param['nu_B']*(LLu*u)
    #bidiff_v = param['nu_B']*(LLv*v)                      
    
    # symmetric stress tensor S = (S11, S12, S12, -S11), store only S11, S12
    S = (Gux*u - Gvy*v,G2vx*v + G2uy*u)
    hS = (h*S[0],h_q*S[1])

    diff_u = (GTx*hS[0] + Gqy*hS[1]) / h_u
    diff_v = (Gqx*hS[1] - GTy*hS[0]) / h_v

    # biharmonic stress tensor R = (R11, R12, R12, -R11), store only R11, R12
    R = (Gux*diff_u - Gvy*diff_v, G2vx*diff_v + G2uy*diff_u)
    hR = (h*R[0],h_q*R[1])
    
    bidiff_u = param['nu_B']*(GTx*hR[0] + Gqy*hR[1]) / h_u
    bidiff_v = param['nu_B']*(Gqx*hR[1] - GTy*hR[0]) / h_v
    
    ## RIGHT-HAND SIDE: ADD TERMS
    rhs_u = adv_u - GTx*p + Fx/h_u - bidiff_u - bfric_u
    rhs_v = adv_v - GTy*p - bidiff_v - bfric_v
    rhs_eta = -(Gux*U + Gvy*V)
    
    return rhs_u, rhs_v, rhs_eta

def ALadvection(q,U,V):
    """ Arakawa and Lamb,1981 advection terms. See interpolation.py for further information. """
    
    AL1q = AL1*q
    AL2q = AL2*q
    
    adv_u = Seur*((ALeur*q)*U) + Seul*((ALeul*q)*U) +\
            Sau*(AL1q[indx_au]*V) + Sbu*(AL2q[indx_bu]*V) +\
            Scu*(AL2q[indx_cu]*V) + Sdu*(AL1q[indx_du]*V)

    adv_v = Spvu*((ALpvu*q)*V) + Spvd*((ALpvd*q)*V) -\
            Sav*(AL1q[indx_av]*U) - Sbv*(AL2q[indx_bv]*U) -\
            Scv*(AL2q[indx_cv]*U) - Sdv*(AL1q[indx_dv]*U)
    
    return adv_u, adv_v