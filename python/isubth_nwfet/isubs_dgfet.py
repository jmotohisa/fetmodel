#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Analysis of Short channel DG MOSFET: subthreshold slope, Vth rolloff

References
[1] Liang, X., & Taur, Y. (2004). A 2-D analytical solution for SCEs in DG MOSFETs.
IEEE Transactions on Electron Devices, 51(9), 1385–1391.
https://doi.org/10.1109/TED.2004.832707

"""

import numpy as np
import math
from scipy import optimize
from scipy.special import jv
from scipy.special import yv
from scipy import integrate
import scipy.constants as const
import matplotlib.pyplot as plt
import fetmodel




def func_for_findroot_sce_dgfet0(invlambda1, tOX0, epsOX, epsS):
    """
    Function to find invlambda1 (Eq.(8) in the Ref.[1])
    Normalized with thickness of channel tSi (tOX0=tOX/tSi)
    """
    lhs = epsS * math.tan(math.pi*tOX0*invlambdad1n)
    rhs = epsOX * math.cot(math.pi/2*invlambda1)
    return(lhs-rhs)


def k1_sce_nwfet0(tOX0, epsOX, epsS):
    """
    solution of Eq.(8) in Ref. [1]: tOX0 is normzlied with radius
    """
    invlambda1 = optimize.root_scalar(
        func_for_findroot_sce_dgfet0, args=(tOX0, epsOX, epsS),
        x0=0.8, x1=2)
    return(1/invlambda1.root)


def k1_sce_dgfet(p):
    """
    solution of Eq.(8) in Ref.[1]
    """
    k1 = k1_sce_nwfet0(p.tox/(p.radius*2), p.eps_ox, p.eps_semi)
    return(k1/(p.radius*2))


def coefsb1_dgfet(k, tSi, tOX, epsOX, epsS):
    """
    coefficient E in Eq.(14) of Ref.[1]
    """
    k1 = k*radius
    k2 = k*(radius+tOX)
    A = jv(0, k1)/((jv(0, k1)*yv(0, k2) - jv(0, k2)*yv(0, k1)))
    B = jv(1, k1)/(jv(1, k1)*yv(0, k2) - jv(0, k2)*yv(1, k1))
    C = (1-epsOX/epsS)*jv(0, k1)**2 + (1-epsS/epsOX)*jv(1, k1)**2
    E = math.pi**2*jv(0, k1)*epsOX/(2*math.log(1+tOX/radius)
                                    * epsS*(A*B+(math.pi/2*k1)**2*C))
    return E


def coefsE(k, p):
    """
    coefficient E in Eq.(14) of Ref.[1]
    """
    return(coefsE0(k, p.radius, p.tox, p.eps_ox, p.eps_semi))


def func_psi(rho, y, Vgs, Vds, p, k1):
    """
    potential psi (Eq. (16) of Ref.[1]
    """
    Eg = p.Eg
    Vbi = Eg/2
    E = coefsE(k1, p)
    b1 = E*(Vbi+p.dphi-Vgs)
    c1 = E*(Vds+Vbi+p.dphi-Vgs)
    psi0 = Vgs - p.dphi + (b1*math.sinh(k1*(p.Lg-y)) +
                           c1*math.sinh(k1*y))/math.sinh(k1*p.Lg)*jv(0, k1*rho)
    return psi0


def func_for_integ_rho_psi(rho, y, Vgs, Vds, p, k1):
    return(math.exp(func_psi(rho, y, Vgs, Vds, p, k1)*const.elementary_charge/(p.temp*const.Boltzmann))*rho)


def func_integd_rho_psi(y, Vgs, Vds, p, k1):
    res = integrate.quad(func_for_integ_rho_psi, 0,
                         p.radius, args=(y, Vgs, Vds, p, k1))
    return(1/res[0])


def func_integd_psi(Vgs, Vds, p, k1):
    res = integrate.quad(func_integd_rho_psi, 0, p.Lg, args=(Vgs, Vds, p, k1))
    return(res[0])


def Vth_rolloff_nwfet0(Vds, Vgs, k1, E, p):
    """
    Vth rolloff (Eqs. (24-26) in Ref. [1])
    """
    Eg = p.Eg
    Vbi = Eg/2
    beta = const.elementary_charge/(p.temp*const.Boltzmann)
    beta1 = 1/beta
    b1 = E*(Vbi+p.dphi-Vgs)
    c1 = E*(Vds+Vbi+p.dphi-Vgs)
    D0 = math.sqrt(b1*c1)*math.exp(-k1*p.Lg/2)
    D1 = math.sqrt(beta*D0)*k1
    yc = p.Lg/2 - math.log(c1/b1)/(2*k1)
    dVth1 = -2*D0
    dVth2 = beta1*math.log(math.sqrt(math.pi)/(2*D1*p.Lg)
                           * (math.erf(D1*(p.Lg-yc))+math.erf(D1*yc)))
    dVth3 = beta1*math.log((D1*p.radius)**2 /
                           (2*(1-math.exp(-(D1*p.radius)**2/2))))
    return dVth1+dVth2+dVth3


def Vth_rolloff_nwfet(Vds, Vgs, p):
    """
    Vth rolloff (Eqs. (24-26) in Ref. [1])
    """
    k1 = k1_sce_nwfet(p)
    E = coefsE(k1, p)
    return Vth_rolloff_nwfet0(Vds, Vgs, k1, E, p)


def Ids_isubsth(Vds, Vgs, p):
    """
    Current Ids (Eq. (19))
    """
    k1 = k1_sce_nwfet(p)
    beta = const.elementary_charge/(p.temp*const.Boltzmann)
    ids = 2*math.pi*p.mue*const.Boltzmann*p.temp*p.ni * \
        (1-math.exp(-Vds*beta))/func_integd_psi(Vgs, Vds, p, k1)
    return(ids)


def SS_subs_nwfet0(Vds, Vgs, k1, E, p):
    """
    Analytic formula for SS (Eq. (31))
    """
    Eg = p.Eg
    Vbi = Eg/2
    dpdv = 1-2*E*(Vbi-p.dphi+Vds/2-Vgs)*math.exp(-k1*p.Lg/2) / \
        (math.sqrt((Vbi+p.dphi-Vgs)*(Vds+Vbi+p.dphi-Vgs)))
    v = p.temp*const.Boltzmann/const.elementary_charge*math.log(10.)
    # print(1/dpdv, v)
    return v/dpdv
    # return 1/dpdv


def SS_subs_nwfet(Vds, Vgs, p):
    """
    Analytic formula for SS (Eq. (31))
    """
    k1 = k1_sce_nwfet(p)
    E = coefsE(k1, p)
    return SS_subs_nwfet0(Vds, Vgs, k1, E, p)

def DIBL_nwfet0(Vds1,Vds2,Vgs,k1,E,p):
    """
    DIBL : delta Vt(Vds1) -delta Vt(Vds2)
    """
    return Vth_rolloff_nwfet0(Vds1, Vgs, k1, E, p) - Vth_rolloff_nwfet0(Vds2, Vgs, k1, E, p)

def DIBL_nwfet(Vds1,Vds2,Vgs,p):
    """
    DIBL : delta Vt(Vds1) -delta Vt(Vds2)
    """
    return Vth_rolloff_nwfet(Vds1, Vgs, p) - Vth_rolloff_nwfet(Vds1, Vgs, p)
    
if __name__ == '__main__':
    p = fetmodel.param_cMOSFET()
    p.radius = 10e-9
    p.Lg = 30e-9
    p.eps_semi = 11.6
    p.Eg = 1.12
    p.Rs = 0
    p.Rd = 0
    p.temp = 300.0
    p.ni = 1.45e16
    p.dphi = 0
    p.tox = 1.5e-9
    p.eps_ox = 3.9
    p.mue = 0.03

    k1 = k1_sce_nwfet(p)
    print("k1=", k1, "lambda (nm)=", math.pi/k1*1e9)
    print(coefsE(k1, p))
    Vds = 1
    Vgs = 0.1
    print(Ids_isubsth(Vds, Vgs, p)*p.Lg*1e9)
    Vgs = 0
    print(SS_subs_nwfet(Vds, Vgs, p))

    p.tox = 3e-9
    p.radius = 5e-9
    k = np.arange(0.5, 5, 0.01)
    fval = np.empty_like(k)
    for i, k0 in enumerate(k):
        fval[i] = func_for_findroot_sce_nwfet0(
            k0, p.tox/p.radius, p.eps_ox, p.eps_semi)

    plt.plot(k, fval)
    plt.show()

    Vds = 1
    Vgs = np.arange(0, 0.8, 0.01)
    Ids1 = np.empty_like(Vgs)
    Ids2 = np.empty_like(Vgs)
    p.Lg = 30e-9
    # for i, Vgs0 in enumerate(Vgs):
    #     Ids1[i] = Ids_isubsth(Vds, Vgs0, p)*p.Lg*1e9

    # p.Lg = 70e-9
    # for i, Vgs0 in enumerate(Vgs):
    #     Ids2[i] = Ids_isubsth(Vds, Vgs0, p)*p.Lg*1e9

    # fig, ax = plt.subplots()
    # ax.plot(Vgs, Ids1, label='Ids1')
    # ax.plot(Vgs, Ids2, label='Ids2')
    # ax.set_yscale("log")
    # plt.xlabel('Gate Voltage (V)')
    # plt.ylabel('Drain Current*L (A nm)')
    # plt.legend(loc='best')
    # plt.show()

    Lg = np.arange(20.0, 100.0, 2.0)
    dVth1 = np.empty_like(Lg)
    dVth2 = np.empty_like(Lg)
    ss1 = np.empty_like(Lg)
    ss2 = np.empty_like(Lg)
    Vds = 1
    Vgs = 0

    p.radius = 5.e-9
    p.tox = 3e-9
    k1 = k1_sce_nwfet(p)
    e1 = coefsE(k1, p)
    print(k1*p.radius, math.pi/k1*1e9, e1)
    for i, Lg0 in enumerate(Lg):
        p.Lg = Lg0*1e-9
        dVth1[i] = Vth_rolloff_nwfet0(Vds, Vgs, k1, e1, p)
        ss1[i] = SS_subs_nwfet0(Vds, Vgs, k1, e1, p)*1e3

    p.radius = 10e-9
    p.tox = 1.5e-9
    k1 = k1_sce_nwfet(p)
    e1 = coefsE(k1, p)
    print(k1*p.radius, math.pi/k1*1e9, e1)
    for i, Lg0 in enumerate(Lg):
        p.Lg = Lg0*1e-9
        dVth2[i] = Vth_rolloff_nwfet0(Vds, Vgs, k1, e1, p)
        ss2[i] = SS_subs_nwfet0(Vds, Vgs, k1, e1, p)*1e3

    plt.plot(Lg, dVth1, label='R=5nm')
    plt.plot(Lg, dVth2, label='R=10nm')
    plt.xlabel('L (m)')
    plt.ylabel('Vth shifft (V)')
    plt.legend(loc='best')
    plt.show()

    plt.plot(Lg, ss1, label='R=5nm')
    plt.plot(Lg, ss2, label='R=10nm')
    plt.xlabel('L (m)')
    plt.ylabel('SS (mV/dec)')
    plt.legend(loc='best')
    plt.show()

    Vds1=0.05
    Vds2=1
    DIBL1 = np.empty_like(Lg)
    # DIBL2 = np.empty_like(Lg)
    for i,Lg0 in enumerate(Lg):
        p.Lg = Lg0*1e-9
        DIBL1[i]= DIBL_nwfet0(Vds1,Vds2,0,k1,e1,p)
        # DIBL2[i]= DIBL_nwfet0(Vds1,Vds2,1,k1,e1,p)
        
    plt.plot(Lg, DIBL1, label='Vgs=0')
    # plt.plot(Lg, DIBL2, label='Vgs=1')
    plt.xlabel('L (m)')
    plt.ylabel('DIBL (V)')
    plt.legend(loc='best')
    plt.show()

    # radius = np.arange(0.010, 0.1, 0.005)
    # l3 = np.empty_like(radius)
    # l2 = np.empty_like(radius)
    # for i, r0 in enumerate(radius):
    #     p.radius = r0*1e-6
    #     l3[i] = lambda3_naturallength(p)
    #     l2[i] = lambda2_naturallength(p)

    # plt.plot(radius, l3, label='lambda3')
    # plt.plot(radius, l2, label='lambda2')
    # plt.legend(loc='best')
    # plt.show()
