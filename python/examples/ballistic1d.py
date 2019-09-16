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
import pdb


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


def check_func_for_E0_rect1dNP(Vds,Vgs,p,left=-0.2,right=0):
    """
    Plot function (Python implementation) to find E0
    Nonparabolic band
    """
    left0 = -(p.alpha_D*Vds+p.alpha_G*Vgs) - 0.2
    left2=min([left,left0])
    right0=-(p.alpha_D*Vds+p.alpha_G*Vgs)
    right2=max([right,right0])
    ene0_list=np.linspace(left2,right2,endpoint=True)
    val=np.empty_like(ene0_list)
    for i,ene0 in enumerate(ene0_list):
        val[i]=func_for_findroot_E0_rect1dNP(ene0, Vds, Vgs, p)

    plt.plot(ene0_list,val)
    plt.show()
    return val


def check_func_for_E0_rect1dNP_fetmodel(Vds,Vgs,p,left=-0.2,right=0):
    """
    Plot function fo find E0 based on fetmodel.density1d_rect1dNP_all0
    Nonparabolic band
    """
    left0 = -(p.alpha_D*Vds+p.alpha_G*Vgs) - 0.2
    left2=min([left,left0])
    right0=-(p.alpha_D*Vds+p.alpha_G*Vgs)
    right2=max([right,right0])
    ene0_list=np.linspace(left2,right2,endpoint=True)
    val=np.empty_like(ene0_list)
    for i,ene0 in enumerate(ene0_list):
        val[i]=fetmodel.func_for_findroot_E0_rect1dNP(ene0, Vds, Vgs, p)

    plt.plot(ene0_list,val)
    plt.show()
    return val


def E0_rect1dNP_root(Vds,Vgs,p,left=-0.2,right=0):
    """
    Get top of the barrier (Python implementation)
    Nonparabolic band
    """
    left0 = -(p.alpha_D*Vds+p.alpha_G*Vgs) - 0.2
    left2=min([left,left0])
    right0=-(p.alpha_D*Vds+p.alpha_G*Vgs)
    right2=max([right,right0])
    e0 = optimize.root_scalar(func_for_findroot_E0_rect1dNP,
                              args=(Vds, Vgs, p), x0=left2, x1=right2)
    if e0.converged==True:
        return e0.root
    else:
        print("EFs convergence error !")
        print(e0)
        return 0
# def get_E0(p, Vds, Vgs):
#     e0 = optimize.root_scalar(func_e0_find, args=(p, Vds, Vgs), x0=-0.1, x1=1)
#     return e0.root


def Ids_ballistic1d_rect1dNP(Vds, Vgs, p, EFs,left=-0.2,right=0):
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


def E0_rect1d_root(Vds,Vgs,p,left=-0.2,right=0):
    """
    Get top of the barrier (Python implementation)
    Parabolic band
    """
    left0 = -(p.alpha_D*Vds+p.alpha_G*Vgs) - 0.2
    left2=min([left,left0])
    right0=-p.alpha_D*Vds-p.alpha_G*Vgs;
    right2=max([right,right0])
    e0 = optimize.root_scalar(func_for_findroot_E0_rect1d,
                              args=(Vds, Vgs, p), x0=left2, x1=right2)
    if e0.converged==True:
        return e0.root
    else:
        print("EFs convergence error !")
        print(e0)
        return 0


def Ids_ballistic1d_rect1d(Vds, Vgs, p, EFs,left=-0.2,right=0):
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
    EFermi=0
    # Eg = 0.36
    # epsOX = 8.5
    # epsS = 8.9
    # ems = 0.2
    # tOX = 20e-9
    Eg = 0.36
    epsOX = 20
    epsS = 15.15
    ems = 0.023
    tOX = 3e-9
    temperature = 300
    W1 = 10e-9
    W2 = 8e-9
    alpha = fetmodel.alpha_NP(Eg, ems)
    Cox = fetmodel.Cox_rect(epsOX, tOX, W1, W2)
    Cc = fetmodel.Cc_rect(epsS, W1, W2)
    # alpha_D = 0
    # alpha_G = 1
    print(Cox,Cc)
    p=fetmodel.parameters_ballistic(EFermi=EFermi,
                                    alpha=alpha,
                                    Ceff=Cox*Cc/(Cox+Cc),
                                    ems=ems,
                                    W1=W1,
                                    W2=W2,
                                    nmax=3,
                                    mmax=4)
    p.output()
    print("Test of density1d: nonparabolic band")

    # p = fetmodel.param_ballistic_new()
    # p.EFermi = 0
    # p.ems = ems
    # p.alpha = alpha
    # p.W1 = W1
    # p.W2 = W2
    # p.alpha_D = alpha_D
    # p.alpha_G = alpha_G
    # p.Ceff = Cox*Cc/(Cox+Cc)
    # p.temp = temperature
    # p.nmax = 3
    # p.mmax = 3

    ### 1DEG density
    EFermi=np.linspace(-0.1,0.5,endpoint=True)
    n1=np.empty_like(EFermi)
    n2=np.empty_like(EFermi)
    n3=np.empty_like(EFermi)

    # ## parabollic band
    # for i,EFermi0 in enumerate(EFermi):
    #     # pdb.set_trace()
    #     p.EFermi=EFermi0
    #     n1[i]=fetmodel.density1d_rect1d_all0(EFermi0, p.ems, p.temp, W1, 8e-9, p.nmax, p.mmax)
    #     n2[i]=fetmodel.density1d_rect1d_all0(EFermi0, p.ems, p.temp, W1, 16e-9, p.nmax, p.mmax)
    #     n3[i]=fetmodel.density1d_rect1d_all0(EFermi0, p.ems, p.temp, W1, 24e-9, p.nmax, p.mmax)

    # fig = plt.figure()
    # plt.plot(EFermi,n1,label="8nm, parabollic")
    # plt.plot(EFermi,n2,label="16nm, parabollic")
    # plt.plot(EFermi,n3,label="24nm, parabollic")
    # plt.xlabel('Fermi Energy (eV)')
    # plt.ylabel('1DEG Density (cm$^{-1}$)')
    # plt.legend(loc='best')
    
    ## nonparabollic band
    for i,EFermi0 in enumerate(EFermi):
        n1[i]=fetmodel.density1d_rect1dNP_all0(
            EFermi0, p.alpha, p.ems, p.temp, W1, 8e-9, p.nmax, p.mmax)
        n2[i]=fetmodel.density1d_rect1dNP_all0(
            EFermi0, p.alpha, p.ems, p.temp, W1, 16e-9, p.nmax, p.mmax)
        n3[i]=fetmodel.density1d_rect1dNP_all0(
            EFermi0, p.alpha, p.ems, p.temp, W1, 24e-9, p.nmax, p.mmax)

    fig = plt.figure()
    plt.plot(EFermi,n1,label="8nm, nonparabollic")
    plt.plot(EFermi,n2,label="16nm, nonparabollic")
    plt.plot(EFermi,n3,label="24nm, nonparabollic")
    plt.legend(loc='best')
    plt.xlabel('Fermi Energy (eV)')
    plt.ylabel('1DEG Density (cm$^{-1}$)')
    plt.show()

    ### EFeremi dependencde of the Top of the barrier
    Vds=0
    Vgs=0.2
    left=-0.3
    right=0.0
    # check_func_for_E0_rect1dNP(Vds,Vgs,p,left=-0.1,right=0.1)
    ene0_list=np.linspace(left,right,endpoint=True)
    val1=np.empty_like(ene0_list)
    val2=np.empty_like(ene0_list)
    val3=np.empty_like(ene0_list)

    p.EFermi=-0.2
    for i,ene0 in enumerate(ene0_list):
        val1[i]=func_for_findroot_E0_rect1dNP(ene0,Vds,Vgs,p)
        # val1[i]=func_for_findroot_E0_rect1d(ene0,Vds,Vgs,p)

    p.EFermi=0
    for i,ene0 in enumerate(ene0_list):
        val2[i]=func_for_findroot_E0_rect1dNP(ene0,Vds,Vgs,p)
        # val2[i]=func_for_findroot_E0_rect1d(ene0,Vds,Vgs,p)

    p.EFermi=0.2
    for i,ene0 in enumerate(ene0_list):
        val3[i]=func_for_findroot_E0_rect1dNP(ene0,Vds,Vgs,p)
        # val3[i]=func_for_findroot_E0_rect1d(ene0,Vds,Vgs,p)
        

    fig = plt.figure()
    plt.plot(ene0_list,val1,label="-0.2")
    plt.plot(ene0_list,val2,label="0")
    plt.plot(ene0_list,val3,label="0.2")
    plt.hlines([0], left,right, "blue", linestyles='dashed')
    plt.legend(loc='best')
    
    Vds = 0.
    Vgs = np.linspace(0, 1,endpoint=True)
    val1 = np.empty_like(Vgs)
    val2 = np.empty_like(Vgs)
    val3 = np.empty_like(Vgs)
    
    p.EFermi=-0.2
    for i, Vgs0 in enumerate(Vgs):
       val1[i] = E0_rect1dNP_root(Vds,Vgs0,p)
        # val1[i] = E0_rect1d_root(Vds,Vgs0,p)
#         val1[i] = fetmodel.E0_rect1dNP_root(Vds,Vgs0,p)
    
    p.EFermi=0
    for i, Vgs0 in enumerate(Vgs):
       val2[i] = E0_rect1dNP_root(Vds,Vgs0,p)
        # val2[i] = E0_rect1d_root(Vds,Vgs0,p)
    
    p.EFermi=0.2
    for i, Vgs0 in enumerate(Vgs):
       val3[i] = E0_rect1dNP_root(Vds,Vgs0,p)
        # val3[i] = E0_rect1d_root(Vds,Vgs0,p)

    fig = plt.figure()
    plt.plot(Vgs, val1,label="-0.2")
    plt.plot(Vgs, val2,label="0")
    plt.plot(Vgs, val3,label="0.2")
    plt.xlabel('Gate Voltage (V)')
    plt.ylabel('Top of the barrier height (eV)')
    plt.legend(loc='best')
    plt.show()

    ## current
    Vds=0.5
    Vgs = np.linspace(0, 1,endpoint=True)
    val1 = np.empty_like(Vgs)
    val2 = np.empty_like(Vgs)
    val3 = np.empty_like(Vgs)
    
    p.EFermi=-0.2
    for i, Vgs0 in enumerate(Vgs):
        val1[i] = Ids_ballistic1d_rect1dNP(Vds, Vgs0, p, 0)/(2*(p.W1+p.W2))
        # val1[i] = Ids_ballistic1d_rect1d(Vds, Vgs0, p, 0)/(2*(p.W1+p.W2))
    p.EFermi=0
    for i, Vgs0 in enumerate(Vgs):
        val2[i] = Ids_ballistic1d_rect1dNP(Vds, Vgs0, p, 0)/(2*(p.W1+p.W2))
        # val2[i] = Ids_ballistic1d_rect1d(Vds, Vgs0, p, 0)/(2*(p.W1+p.W2))
    p.EFermi=0.2
    for i, Vgs0 in enumerate(Vgs):
        val3[i] = Ids_ballistic1d_rect1dNP(Vds, Vgs0, p, 0)/(2*(p.W1+p.W2))
        # val3[i] = Ids_ballistic1d_rect1d(Vds, Vgs0, p, 0)/(2*(p.W1+p.W2))

    fig, ax = plt.subplots()
    ax.plot(Vgs, val1, label='EF=-0.2 eV')
    ax.plot(Vgs, val2, label='EF=0 eV')
    ax.plot(Vgs, val3, label='EF=0.2 eV')
    ax.set_yscale("log")
    plt.hlines([100e-3], min(Vgs),max(Vgs), "blue", linestyles='dashed')
    plt.xlabel('Gate Voltage (V)')
    plt.ylabel('Drain Current (uA/um)')
    plt.legend(loc='best')
        
    Vds=0.5
    p.EFermi=0
    Vds = np.arange(0, 1, 0.01)
    Ids1 = np.empty_like(Vds)
    Ids2 = np.empty_like(Vds)
    Vgs=0.5
    for i, Vds0 in enumerate(Vds):
        Ids1[i] = Ids_ballistic1d_rect1dNP(Vds0, Vgs, p, 0)/(2*(p.W1+p.W2))
        # Ids1[i] = Ids_ballistic1d_rect1d(Vds0, Vgs, p, 0)/(2*(p.W1+p.W2))
    Vgs=1
    for i, Vds0 in enumerate(Vds):
        Ids2[i] = Ids_ballistic1d_rect1dNP(Vds0, Vgs, p, 0)/(2*(p.W1+p.W2))
        # Ids2[i] = fetmodel.Ids_ballistic1d_rect1dNP(Vds0, Vgs, p, 0)
        # Ids2[i] = Ids_ballistic1d_rect1d(Vds0, Vgs, p, 0)/(2*(p.W1+p.W2))

    fig = plt.figure()
    plt.plot(Vds, Ids1, label='Vgs=0.5 V')
    plt.plot(Vds, Ids2, label='Vgs=1 V')
    plt.xlabel('Drain Voltage (V)')
    plt.ylabel('Drain Current (uA/um)')
    plt.legend(loc='best')

    plt.show()
