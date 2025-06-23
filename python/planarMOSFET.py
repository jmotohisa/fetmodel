#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Python module for planar MOSFET

import numpy as np
import scipy.constants as const
import fetmodel
from scipy import optimize


def psiB_func(p):
    """
    p class of parameteres_plMOSFET
    return Differece of Fermi Energy and intrinsic Fermi Energy in p-type semiconductor
    phi_F - phi_Fi
    """
    return (p.temp*const.k/const.e)*np.log(p.NA/p.ni)


def Vth_body(Vbs, p):
    """
    Vbs source to substrate voltage
    p class of parameteres_plMOSFET
    return Threshold voltage with body effect
    """
    psiB = psiB_func(p)
    return 2*psiB+np.sqrt(2*const.epsilon_0*p.eps_semi*const.e*p.NA*(2*psiB-Vbs))/p.Cox


def Vth_plMOSFET(p):
    """
    Vbs source to substrate voltage
    p class of parameteres_plMOSFET
    return Threshold voltage
    """
    psiB = psiB_func(p)
    return 2*psiB+np.sqrt(2*const.epsilon_0*p.eps_semi*const.e*p.NA*2*psiB)/p.Cox


def Vdssat1_plMOSFET(Vgs, p):
    """
    return saturation drain voltage
    """
    psiB = psiB_func(p)
    c = (p.eps_semi*const.epsilon_0*const.e*p.NA)/(p.Cox**2)
    return Vgs-2*psiB+c*(1-np.sqrt(1+2*(Vgs)/c))


def Ids0_plMOSFET(Vgs, Vds, Vth, p):
    """
    Drain current of MOSFET: analitic model (The most basic model)
    Vgs Gate volatage
    Vds Drain voltage
    Vth Threshold Voltage
    p  class of parameters_plMOSFET
    """
    Vdssat = Vgs-Vth
    c = p.mue*p.Cox/p.Lg
    if (Vgs > Vth):
        if (Vds < Vdssat):
            return c*((Vgs-Vth)*Vds-Vds**2/2)
        else:
            return c/2*(Vgs-Vth)**2
    else:
        return 0


def Ids1_func0(Vgs, Vds, p):
    """
    Drain current of MOSFET
    Vgs Gate volatage
    Vds Drain voltage
    p  class of parameters_plMOSFET
    """
    psiB = psiB_func(p)
    c = p.mue*p.Cox/p.Lg
    return c*((Vgs-2*psiB)*Vds-Vds**2/2-2/3*np.sqrt(2*p.eps_semi*const.epsilon_0*const.e*p.NA)/p.Cox*((2*psiB+Vds)**(1.5)-(2*psiB)**(1.5)))


def Ids1_plMOSFET(Vgs, Vds, p):
    """
    Drain current of MOSFET
    Vgs Gate volatage
    Vds Drain voltage
    p  class of parameters_plMOSFET
    """
    Vdssat = Vdssat1_plMOSFET(Vgs, p)
    Vth0 = Vth_plMOSFET(p)
    if (Vgs > Vth0):
        if (Vds < Vdssat):
            return Ids1_func0(Vgs, Vds, p)
        else:
            return Ids1_func0(Vgs, Vdssat, p)
    else:
        return 0


def Ids2_plMOSFET(Vgs, Vds, p):
    """
    Drain current of MOSFET in Regional approximation
    Vgs Gate volatage
    Vds Drain voltage
    p  class of parameters_plMOSFET
    """
    Vth0 = Vth_plMOSFET(p)
    psiB = psiB_func(p)
    m = 1+np.sqrt(p.eps_semi*const.epsilon_0*const.e*p.NA/(4*psiB))/p.Cox
    Vdssat = (Vgs-Vth0)/m
    if (Vgs > Vth0):
        if (Vds < Vdssat):
            return Ids2_func0(Vgs, Vds, p)
        else:
            return Ids2_func0(Vgs, Vdssat, p)
    else:
        return 0


def Ids2_func0(Vgs, Vds, p):
    """
    Drain current of MOSFET in Regional approximation in the triode region
    Vgs Gate volatage
    Vds Drain voltage
    p  class of parameters_plMOSFET
    """
    Vth0 = Vth_plMOSFET(p)
    psiB = psiB_func(p)
    m = 1+np.sqrt(p.eps_semi*const.epsilon_0*const.e*p.NA/(4*psiB))/p.Cox
    c = p.mue*p.Cox/p.Lg
    return c*((Vgs-Vth0)*Vds-Vds**2*m/2)


def func_psiS_findroot_plMOSFET(psiS, p, Vgs, Vds):
    """
    """
    kBT = const.k*p.temp
    a = np.sqrt(2*p.eps_semi*const.epsilon_0*kBT*p.NA)/p.Cox
    b = const.e*psiS/kBT+p.ni**2/(p.NA**2)*np.exp(const.e*(psiS-Vds)/kBT)
    return psiS+a*np.sqrt(b)-Vgs


def func_psiS_findroot_prime_plMOSFET(psiS, p, Vgs, Vds):
    kBT = const.k*p.temp
    a = np.sqrt(2*p.eps_semi*const.epsilon_0*kBT*p.NA)/p.Cox
    b = const.e*psiS/kBT+p.ni**2/(p.NA**2)*np.exp(const.e*(psiS-Vds)/kBT)
    return 1+a/np.sqrt(b)*const.e/kBT*(1+p.ni**2/(p.NA**2)*np.exp(const.e*(psiS-Vds)/kBT))


def func_psiS_findroot1_plMOSFET(psiS, p, Vgs):
    kBT = const.k*p.temp
    psiB = psiB_func(p)
    a = np.sqrt(2*p.eps_semi*const.epsilon_0*const.e*p.NA*psiS)/p.Cox
    return psiB+a-Vgs


def find_psiS_plMOSFET(p, Vgs, Vds):
    x0 = psiB_func(p)
    e0 = optimize.root_scalar(func_psiS_findroot_plMOSFET, method='newton',
                              fprime=func_psiS_findroot_prime_plMOSFET,
                              args=(p, Vgs, Vds), x0=x0)
    return e0.root
