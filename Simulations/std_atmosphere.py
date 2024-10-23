# Standard atmosphere for Troposphere layer
def geometric_altitud(H, r=6356766):
    h = r*H/(r-H)
    return h

def temperature(H):
    if H <= 11000: 
        Tb = 288.15
        Hb = 0
        beta = -0.0065
    else:
        Tb = 216.65
        Hb = 11000
        beta = 0
    T = Tb + beta*(H - Hb)
    return T

def pressure(H):
    import numpy as np
    R = 287.05287
    g0=9.80665
    if H <= 11000:
        pb = 101325
        Hb = 0
        Tb = 288.15
        beta = -0.0065
        p = pb*(1 + beta*(H - Hb)/Tb)**(-g0/(beta*R))
    else:
        pb = 22632
        Hb = 11000
        Tb = 216.15
        p = pb*np.exp(-g0*(H-Hb)/(R*Tb))
    return p

def density(p, T):
    R = 287.05287
    rho = p/(R*T)
    return rho

def sound_speed(T):
    a = 20.046796*T**0.5
    return a

def dynamic_viscosity(T):
    beta_S = 1.458e-6
    S = 110.4
    mu = beta_S*T**(3/2)/(T + S)
    return mu

def cinematic_viscosity(mu, rho):
    nu = mu/rho
    return nu

def specific_heat_cp(T):
    cp = 1061.332 - 0.432819*T + 1.02344e-3*T**2 - 6.47474e-7*T**3 + 1.3864e-10*T**4
    return cp

def thermal_conductivity(T):
    kappa = 2.648151e-3*T**1.5/(T+245.4*10**(-12/T))
    return kappa

def Prandtl_number(mu, cp, kappa):
    Pr = mu*cp/kappa
    return Pr

def simul_conditions(Re, Ma, c):
    Hn_1 = 0
    Hn = 11000
    Eobj = 1e-4
    E = abs(Hn - Hn_1)/abs(Hn)
    while (E>Eobj):
        Tn_1 = temperature(Hn_1)
        pn_1 = pressure(Hn_1)
        rhon_1 = density(pn_1, Tn_1)
        an_1 = sound_speed(Tn_1)
        mun_1 = dynamic_viscosity(Tn_1)
        
        fn_1 = an_1*rhon_1/mun_1 - Re/(Ma*c)
        
        Tn = temperature(Hn)
        pn = pressure(Hn)
        rhon = density(pn, Tn)
        an = sound_speed(Tn)
        mun = dynamic_viscosity(Tn)
        
        fn = an*rhon/mun - Re/(Ma*c)

        Hnm1 = Hn - (Hn - Hn_1)*fn/(fn - fn_1)
        
        E = abs(Hnm1 - Hn)/abs(Hnm1)
        
        Hn_1 = Hn
        Hn = Hnm1
    return Hnm1