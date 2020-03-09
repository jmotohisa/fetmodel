#!/usr/bin/env python
# -*- coding: utf-8 -*-

# sample Python script of cyrlindircal MOSFET
# B. Iniguez et al., IEEE Trans. Elec. Dev. 52, No.8, Auguust 2015 (p.1868)

import fetmodel
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import scipy.constants as const
from scipy import integrate
# import h5py

def qfunc_cMOSFET(qq,V,Vgs,p):
    if(qq<0):
        qq=np.abs(qq)/10.

    Vth = const.Boltzmann*p.temp/const.elementary_charge
    Q0= 4*p.eps_semi*const.epsilon_0/p.radius*Vth
    delta = const.elementary_charge**2*p.ni/(const.Boltzmann*p.temp*p.eps_semi*const.epsilon_0)
    qqq1 = Vgs-p.dphi-V-Vth*math.log(8/(delta*p.radius**2))
    qqq2 =(qq/p.Cox + Vth*(math.log(qq/Q0)+math.log(1+(qq/Q0))));
    return qqq1-qqq2

def qroot_brent(V,Vgs,p):
    e0 = optimize.root_scalar(qfunc_cMOSFET,
                              args=(V, Vgs, p), x0=1e-9,x1=1e-3)
    if e0.converged==True:
        return e0.root
    else:
        print("Q-cMOSFET convergence error !")
        print(e0)
        return 0

def Voff_from_Ioff_cMOSFET(p, Vds, Ids_cutoff, ps, left=-0.5, right=0.5):
    e0 = optimize.root_scalar(
        func_det_Voff, args=(p, Vds, Ids_cutoff,ps), x0=left, x1=right,method='brentq',bracket=[left,right])
    if(e0.converged==True):
        return e0.root
    else:
        print("Vgs convergence error !")
        print(e0)
        return 0

def func_det_Voff(Vgs, p, Vds, Ids_cutoff, ps):
    return Ids_cutoff - fetmodel.func_Ids_cMOSFET(Vds, Vgs, p,ps)/(2*math.pi*p.radius)

def Voff_from_Ioff_cMOSFET2(p, Vds, Ids_cutoff, left=-0.5, right=0.5):
    e0 = optimize.root_scalar(
        func_det_Voff2, args=(p, Vds, Ids_cutoff), x0=left, x1=right,method='brentq',bracket=[left,right])
    if(e0.converged==True):
        return e0.root
    else:
        print("dphi convergence error !")
        print(e0)
        return 0

def func_det_Voff2(Vgs, p, Vds, Ids_cutoff):
    return Ids_cutoff - fetmodel.func_Ids2_cMOSFET(Vds, Vgs, p)/(2*math.pi*p.radius)


if __name__ == '__main__':
    p = fetmodel.param_cMOSFET()
    p.radius = 6.25e-9
    p.Lg = 1e-6
    p.eps_semi = 11.6
    p.Rs = 0
    p.Rd = 0
    p.temp = 300
    p.ni = 1.45e16
    p.dphi = 0
    p.tox = 1.5e-9
    p.eps_ox = 3.9
    p.mue = 0.04
    p.Cox = p.eps_ox*8.85e-12/(p.radius*math.log(1+p.tox/p.radius))
    # p.Cox = fetmodel.Cox_radial(p.eps_ox,p.txo,p.radius)

    qq=1e-3
    V=0.5
    Vgs=0
    print(qfunc_cMOSFET(qq,V,Vgs,p))
    print(qroot_brent(V,Vgs,p))
    print(fetmodel.func_rootfind_Q_cMOSFET(qq,V,Vgs,p))

    Vgs = np.arange(0, 2, 0.01)

    # Fig.2
    Q1 = np.empty_like(Vgs)
    Q2 = np.empty_like(Vgs)
    fetmodel.Qapprox_cMOS(Vgs, Q1, p)
    fetmodel.Q_cMOS(Vgs, Q2, p)
    # Q1 = fetmodel.Qapprox_cMOS(Vgs, p)
    # Q2 = fetmodel.Q_cMOS(Vgs, p)
    
    plt.plot(Vgs, Q1)
    plt.plot(Vgs, Q2)
    plt.show()
    
    # Fig.3
    Ids01 = np.empty_like(Vgs)
    Ids02 = np.empty_like(Vgs)
    fetmodel.Ids_cMOS(Vgs, Ids01, 0.1, p)
    fetmodel.Ids0_cMOS(Vgs, Ids02, 0.1, p)
    plt.plot(Vgs, Ids01)
    plt.plot(Vgs, Ids02)
    plt.show()
    
    # Fig.4
    fetmodel.Ids_cMOS(Vgs, Ids01, 1, p)
    fetmodel.Ids0_cMOS(Vgs, Ids02, 1, p)
    plt.plot(Vgs, Ids01)
    plt.plot(Vgs, Ids02)
    plt.show()

