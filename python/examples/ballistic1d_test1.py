#!/usr/bin/env python
# -*- coding: utf-8 -*-

# test of ballistic FET (using lower level interface)

import fetmodel
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from ballistic1d import *


def e0Vgs(p,Vds,Vgs_list):
    e0=np.empty_like(Vgs_list)
    for i,Vgs0 in enumerate(Vgs_list):
        e0[i]=E0_rect1dNP_root(Vds,Vgs0,p)

    return e0

def set_p(p,Eg,W2,nmax,mmax):
    p.ems = func_ems(Eg)
    p.alpha = fetmodel.alpha_NP(Eg, p.ems)
    p.W2 = W2
    Cox = fetmodel.Cox_rect(epsOX, tOX, p.W1, W2)
    Cc = fetmodel.Cc_rect(epsS, p.W1, W2)
    p.Ceff = Cox*Cc/(Cox+Cc)
    p.nmax = nmax
    p.mmax = mmax

def func_ems(Eg,q20=20):
    return 1/(1+q20/Eg)


def Ids_ballistic1d_rect1dNP_E0(E0, Vds, Vgs, p, EFs):
    nlist = np.arange(1, p.nmax+1, dtype=np.int64)
    mlist = np.arange(1, p.mmax+1, dtype=np.int64)
    cur = 0
    for n in nlist:
        for m in mlist:
            Enmp = fetmodel.Ep_nm_rect1d(p.ems, p.W1, p.W2, int(n), int(m))
            gamma_nm = fetmodel.gamma_nm_NP(Enmp, p.alpha)
            Enm = fetmodel.E_nm_NP(p.alpha, gamma_nm)
            cur1 = func_FD0(EFs-Enm-E0, p.temp)
            cur2 = func_FD0(EFs-Enm-E0-Vds, p.temp)
            cur += cur1-cur2

    return (cur*2*const.elementary_charge/const.h*p.temp*const.Boltzmann)

def IdsVgsE0(Vds,Vgs_list,E0_list,EFs):
    Ids=np.empty_like(Vgs_list)
    for i, Vgs0 in enumerate(Vgs_list):
        Ids[i]=Ids_ballistic1d_rect1dNP_E0(E0_list[i], Vds, Vgs0, p, EFs)

    return Ids


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

    Vds=0.5
    Vgs_list=np.linspace(0.0,0.5,endpoint=True)
    Eg=3.4
    W2=8e-9
    nmax=3
    mmax=3
    set_p(p,Eg,W2,nmax,nmax)
    E01=e0Vgs(p,Vds,Vgs_list)
    Ids1=IdsVgsE0(Vds,Vgs_list,E01,0)/(2*(p.W1+p.W2))

    W2=16e-9
    set_p(p,Eg,W2,nmax,nmax)
    E02=e0Vgs(p,Vds,Vgs_list)
    Ids2=IdsVgsE0(Vds,Vgs_list,E02,0)/(2*(p.W1+p.W2))

    W2=24e-9
    set_p(p,Eg,W2,nmax,nmax)
    E03=e0Vgs(p,Vds,Vgs_list)
    Ids3=IdsVgsE0(Vds,Vgs_list,E02,0)/(2*(p.W1+p.W2))

    fig=plt.figure()
    plt.plot(Vgs_list,E01,label = 'W2=8 nm')
    plt.plot(Vgs_list,E02,label = 'W2=16 nm')
    plt.plot(Vgs_list,E03,label = 'W2=24 nm')
    plt.xlabel('Gate Voltage (V)',fontsize=18)
    plt.ylabel('Top of the barrier height (eV)',fontsize=18)
    plt.tick_params(labelsize=18)
    plt.legend(loc='best',fontsize=18)
    plt.tight_layout()

    fig,ax = plt.subplots()
    ax.plot(Vgs_list,Ids1,label = 'W2=8 nm')
    ax.plot(Vgs_list,Ids2,label = 'W2=16 nm')
    ax.plot(Vgs_list,Ids3,label = 'W2=16 nm')
    ax.set_yscale("log")
    plt.xlabel('Gate Voltage (V)',fontsize=18)
    plt.ylabel('Current (uA/um)',fontsize=18)
    plt.tick_params(labelsize=18)
    plt.legend(loc='best',fontsize=18)
    plt.tight_layout()
    
    plt.show()

    
    # print(fetmodel.E0_rect1d_root(p))
    # print(p.alpha, p.ems, p.temp, p.W1, p.W2, p.nmax, p.mmax, p.Ceff)
    # print(fetmodel.density1d_rect1dNP_all0(
    #     0.1, p.alpha, p.ems, p.temp, p.W1, p.W2, p.nmax, p.mmax))
    
    
    # print(alpha, gamma_nm, Enm, alpha_nm, ems_nm)
    
    Vds = 0
    Vgs = np.arange(0, 1, 0.01)
    e0 = np.empty_like(Vgs)
    
    Vds = np.arange(0, 1, 0.01)
    Ids1 = np.empty_like(Vds)
    Ids2 = np.empty_like(Vds)
        
    # for i, Vds0 in enumerate(Vds):
    #     Ids1[i] = Ids_ballistic1d_rect1dNP(Vds, Vgs_list, p, EFs)/(2*(p.W1+p.W2))
    #     # Ids2[i] = Ids_ballistic1d_rect1dNP(Vds, Vgs_list, p, EFs)/(2*(p.W1+p.W2))
        
    # plt.plot(Vds, Ids1, label='Ids1')
    # # plt.plot(Vds, Ids2, label='Ids2')
    # plt.xlabel('Drain Voltage (V)')
    # plt.ylabel('Drain Current (A)')
    # plt.legend(loc='best')
    # plt.show()
        
