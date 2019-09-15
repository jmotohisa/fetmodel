#!/usr/bin/env python
# -*- coding: utf-8 -*-

# test of ballistic FET (using lower level interface)

import fetmodel
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import scipy.constants as const
from scipy import integrate


def func_for_findroot_E0_rect1dNP(ene0, Vds, Vgs, p):
    """
    Python implementation of the function to find top of the barrier
    based on fetmodel.density1d_rect1dNP_all0
    Nonparabolic band
    """
    n1d_S = fetmodel.density1d_rect1dNP_all0(
        p.EFermi - ene0, p.alpha, p.ems, p.temp, p.W1, p.W2, p.nmax, p.mmax)
    n1d_D = fetmodel.density1d_rect1dNP_all0(
        p.EFermi - ene0 - Vds, p.alpha, p.ems, p.temp, p.W1, p.W2, p.nmax, p.mmax)
    q0 = const.elementary_charge * (n1d_S + n1d_D) / (2 * p.Ceff)
    return ene0 + (p.alpha_D * Vds + p.alpha_G * Vgs - q0)
# def func_e0_find(E0, p, Vds, Vgs):
#     n1d_S = fetmodel.density1d_rect1dNP_all0(
#         p.EFermi - E0, p.alpha, p.ems, p.temp, p.W1, p.W2, p.nmax, p.mmax)
#     n1d_D = fetmodel.density1d_rect1dNP_all0(
#         p.EFermi - E0 - Vds, p.alpha, p.ems, p.temp, p.W1, p.W2, p.nmax, p.mmax)
#     q0 = 1.6e-19 * (n1d_S + n1d_D) / (2 * p.Ceff)
#     return E0 + (p.alpha_D * Vds + p.alpha_G * Vgs - q0)


def check_func_for_E0_rect1dNP(E0start,E0stop,Vds,Vgs,p):
    """
    Plot function (Python implementation) to find E0
    Nonparabolic band
    """
    ene0_list=np.linspace(E0start,E0stop,endpoint=True)
    val=np.empty_like(ene0_list)
    for i,ene0 in enumerate(ene0_list):
        val[i]=func_for_findroot_E0_rect1dNP(ene0, Vds, Vgs, p)

    plt.plot(ene0_list,val)
    plt.show()
    return val


def check_func_for_E0_rect1dNP_fetmodel(E0start,E0stop,Vds,Vgs,p):
    """
    Plot function fo find E0 based on fetmodel.density1d_rect1dNP_all0
    Nonparabolic band
    """
    ene0_list=np.linspace(E0start,E0stop,endpoint=True)
    val=np.empty_like(ene0_list)
    for i,ene0 in enumerate(ene0_list):
        val[i]=fetmodel.func_for_findroot_E0_rect1dNP(ene0, Vds, Vgs, p)

    plt.plot(ene0_list,val)
    plt.show()
    return val


def E0_rect1dNP_root(Vds,Vgs,p,left=-0.1,right=0.1):
    """
    Get top of the barrier (Python implementation)
    Nonparabolic band
    """
    e0 = optimize.root_scalar(func_for_findroot_E0_rect1dNP,
                              args=(p,Vds, Vgs), x0=left, x1=right)
    if e0.converged==True:
        return e0.root
    else:
        print("EFs convergence error !")
        print(e0)
        return 0
# def get_E0(p, Vds, Vgs):
#     e0 = optimize.root_scalar(func_e0_find, args=(p, Vds, Vgs), x0=-0.1, x1=1)
#     return e0.root


def Ids_ballistic1d_rect1dNP(Vds, Vgs, p, EFs,left=-0.1,right=0.1):
    e0=E0_rect1d_root(Vds,Vgs,p,left,right)
    # e00 = optimize.root_scalar(func_e0_find, args=(p, Vgs, Vgs), x0=-0.1, x1=1)
    # e00=get_E0(p, Vgs, Vds)
    # e0=e00.root
    nlist = np.arange(1, p.nmax+1, dtype=np.int64)
    mlist = np.arange(1, p.mmax+1, dtype=np.int64)
    cur = 0
    for n in nlist:
        for m in mlist:
            Enmp = fetmodel.Ep_nm_rect1d(p.ems, p.W1, p.W2, int(n), int(m))
            gamma_nm = fetmodel.gamma_nm_NP(Enmp, p.alpha)
            Enm = fetmodel.E_nm_NP(p.alpha, gamma_nm)
            cur1 = func_FD0(EFs-Enm-e0, p.temp)
            cur2 = func_FD0(EFs-Enm-e0-Vds, p.temp)
            cur += cur1-cur2

    return (cur*2*const.elementary_charge/const.h*p.temp*const.Boltzmann)

####### functions
def func_FD0(ene, temp):
    """
    Fermi-Dirac integral in the zero-th order
    """
    ene0=ene*const.elementary_charge/(const.Boltzmann*temp)
    if(ene0>100):
        return ene0
    else:
        return math.log(1+math.exp(ene0))

def Ep_nm_rect1d(ems, W1, W2, n, m):
    return (math.pi*const.hbar)**2/2*((n/W1)**2+(m/W2)**2)/(ems*const.electron_mass*const.elementary_charge)

def gamma_nm_NP(Enmp, alpha):
    return math.sqrt(1+4*alpha*Enmp)

def E_nm_NP(alpha,gamma_nm):
    return (gamma_nm-1)/(2*alpha)

def ems_nm_NP(ems,gamma_nm):
    return ems*gamma_nm

def alpha_nm_NP(alpha,gamma_nm):
    return alpha/gamma_nm

def dos1D_NP_norm(eta,alpha0):
    return (1+2*alpha0*eta)/math.sqrt(eta*(1+alpha0*eta))

def func_for_integdensity(eta,etaF,alpha0):
    return 1/(1+math.exp(eta-etaF))*dos1D_NP_norm(eta,alpha0)

def density1d_rect1dNP0(EFermi,alpha,ems,temp,W1,W2,n,m):
    Enm = Ep_nm_rect1d(ems, W1, W2, n, m)
    gamma_nm = gamma_nm_NP(Enm, alpha)
    ems_nm = ems_nm_NP(ems,gamma_nm)
    alpha_nm = alpha_nm_NP(alpha,gamma_nm)
    EtaF=(EFermi - Enm)*const.elementary_charge/(temp*const.Boltzmann)
    alpha0 = alpha_nm/const.elementary_charge*(temp*const.Boltzmann)
    res = integrate.quad(func_for_integdensity, 0.001,
                         np.inf, args=(EtaF,alpha0))
    return res[0]*math.sqrt(2*ems*const.electron_mass*temp*const.Boltzmann)/(math.pi*const.hbar)

    
def density1d_rect1dNP_all0(EFermi, alpha, ems, temp, W1, W2, nmax, mmax):
    n0=0.
    for n in np.arange(0,nmax):
        for m in np.arange(0,mmax):
            dens=density1d_rect1dNP0(EFermi,alpha,ems,temp,W1,W2,int(n)+1,int(m)+1)
            print(n,m,dens)
            n0+=dens

    return(n0)

def density1d_rect1dNP_all(p):
    return density1d_rect1dNP_all0(p.EFermi, p.alpha, p.ems, p.temp, p.W1, p.W2, p.nmax, p.mmax)

#####
# parabolic band, rectangular NW
#####

def func_for_findroot_E0_rect1d(ene0, Vds, Vgs, p):
    """
    Python implementation of the function to find top of the barrier
    based on fetmodel.density1d_rect1d_all0
    Parabolic band
    """
    n1d_S = fetmodel.density1d_rect1d_all0(
        p.EFermi - ene0, p.ems, p.temp, p.W1, p.W2, p.nmax, p.mmax)
    n1d_D = fetmodel.density1d_rect1d_all0(
        p.EFermi - ene0 - Vds, p.ems, p.temp, p.W1, p.W2, p.nmax, p.mmax)
    q0 = const.elementary_charge * (n1d_S + n1d_D) / (2 * p.Ceff)
    return ene0 + (p.alpha_D * Vds + p.alpha_G * Vgs - q0)


def check_func_for_E0_rect1d(E0start,E0stop,Vds,Vgs,p):
    """
    Plot function (Python implementation) to find E0
    Parabolic band
    """
    ene0_list=np.linspaace(E0start,E0stop,endpoint=True)
    val=np.empty_like(ene0_list)
    for i,ene0 in enumerate(ene0_list):
        val[i]=func_for_findroot_E0_rect1d(ene0, Vds, Vgs, p)

    plt.plot(ene0_list,val)
    return val


def E0_rect1d_root(Vds,Vgs,p,left=-0.1,right=0.1):
    """
    Get top of the barrier (Python implementation)
    Parabolic band
    """
    e0 = optimize.root_scalar(func_for_findroot_E0_rect1d,
                              args=(Vds, Vgs, p), x0=left, x1=right)
    if e0.converged==True:
        return e0.root
    else:
        print("EFs convergence error !")
        print(e0)
        return 0


def Ids_ballistic1d_rect1d(Vds, Vgs, p, EFs,left=-0.1,right=0.1):
    e0=E0_rect1d_root(Vds,Vgs,p,left,right)
    # e0 = optimize.root_scalar(func_e0_find, args=(p, Vgs, Vgs), x0=-0.1, x1=1)
    # e0=get_E0(p, Vgs, Vds)
    nlist = np.arange(1, p.nmax+1, dtype=np.int64)
    mlist = np.arange(1, p.mmax+1, dtype=np.int64)
    cur = 0
    for n in nlist:
        for m in mlist:
            Enmp = fetmodel.Ep_nm_rect1d(p.ems, p.W1, p.W2, int(n), int(m))
            cur1 = func_FD0(EFs-Enmp-e0, p.temp)
            cur2 = func_FD0(EFs-Enmp-e0-Vds, p.temp)
            cur += cur1-cur2

    return (cur*2*const.elementary_charge/const.h*p.temp*const.Boltzmann)


if __name__ == '__main__':
    Eg = 0.36
    epsOX = 8.5
    epsS = 8.9
    tOX = 20e-9
    temperature = 300
    ems = 0.2
    W1 = 5e-9
    W2 = 8e-9
    alpha = fetmodel.alpha_NP(Eg, ems)
    Cox = fetmodel.Cox_rect(epsOX, tOX, W1, W2)
    Cc = fetmodel.Cc_rect(epsS, W1, W2)
    alpha_D = 0
    alpha_G = 1

    # p = fetmodel.param_ballistic_new()
                                       
    p.ems = ems
    p.alpha = alpha
    p.W1 = W1
    p.W2 = W2

    p.EFermi = 0
    p.alpha_D = alpha_D
    p.alpha_G = alpha_G
    p.Ceff = Cox*Cc/(Cox+Cc)
    p.temp = temperature
    p.nmax = 3
    p.mmax = 3

    EFermi=np.arange(-0.1,0.5,0.01)
    n1=np.empty_like(EFermi)
    n2=np.empty_like(EFermi)
    n3=np.empty_like(EFermi)
    for i,EFermi0 in enumerate(EFermi):
        n1[i]=fetmodel.density1d_rect1dNP_all0(
            EFermi0, p.alpha, p.ems, p.temp, 5e-9, 8e-9, p.nmax, p.mmax)
        n2[i]=fetmodel.density1d_rect1dNP_all0(
            EFermi0, p.alpha, p.ems, p.temp, 5e-9, 16e-9, p.nmax, p.mmax)
        n3[i]=fetmodel.density1d_rect1dNP_all0(
            EFermi0, p.alpha, p.ems, p.temp, 5e-9, 24e-9, p.nmax, p.mmax)

    plt.plot(EFermi,n1,label="8nm")
    plt.plot(EFermi,n2,label="16nm")
    plt.plot(EFermi,n3,label="24nm")
    plt.legend(loc='best')
    plt.show()
    # print(fetmodel.E0_rect1d_root(p))
    # print(p.alpha, p.ems, p.temp, p.W1, p.W2, p.nmax, p.mmax, p.Ceff)
    # print(fetmodel.density1d_rect1dNP_all0(
    #     0.1, p.alpha, p.ems, p.temp, p.W1, p.W2, p.nmax, p.mmax))

    # print(alpha, gamma_nm, Enm, alpha_nm, ems_nm)

    Vgs = 0
    Vds = 0
    e = np.arange(-0.1, 0.1, 0.005)
    e0 = np.empty_like(e)
    for i, e00 in enumerate(e):
        e0[i] = func_e0_find(e00, p, Vgs, Vds)

    plt.plot(e, e0)
    plt.xlabel('Energy (eV)')
    plt.ylabel('func_E0_find')
    plt.show()

    Vds = 0
    Vgs = np.arange(0, 1, 0.01)
    e0 = np.empty_like(Vgs)
    
    for i, Vgs0 in enumerate(Vgs):
        e0[i] = get_E0(p, Vgs0, Vds)
        
    plt.plot(Vgs, e0)
    plt.xlabel('Gate Voltage (V)')
    plt.ylabel('Top of the barrier height (eV)')
    plt.show()


    Vds = np.arange(0, 1, 0.01)
    Ids1 = np.empty_like(Vds)
    Ids2 = np.empty_like(Vds)
    
    for i, Vds0 in enumerate(Vds):
        Ids1[i] = func_current1D(-0.1, Vds0, p, 0)
        Ids2[i] = func_current1D(0, Vds0, p, 0)

    plt.plot(Vds, Ids1, label='Ids1')
    plt.plot(Vds, Ids2, label='Ids2')
    plt.xlabel('Drain Voltage (V)')
    plt.ylabel('Drain Current (A)')
    plt.legend(loc='best')
    plt.show()
