#!/usr/bin/env python
# -*- coding: utf-8 -*-

import fetmodel
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import scipy.constants as const
from ballistic1d import *


if __name__ == '__main__':
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
    p=fetmodel.parameters_ballistic(alpha=alpha,
                                    Ceff=Cox*Cc/(Cox+Cc),
                                    ems=ems,
                                    W1=W1,
                                    W2=W2,
                                    nmax=3,
                                    mmax=4)
    p.output()
    print("Test of density1d: parabolic band")
    # p = fetmodel.param_ballistic()
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

    ## parabollic band
    for i,EFermi0 in enumerate(EFermi):
        # pdb.set_trace()
        p.EFermi=EFermi0
        n1[i]=fetmodel.density1d_rect1d_all0(EFermi0, p.ems, p.temp, W1, 8e-9, p.nmax, p.mmax)
        n2[i]=fetmodel.density1d_rect1d_all0(EFermi0, p.ems, p.temp, W1, 16e-9, p.nmax, p.mmax)
        n3[i]=fetmodel.density1d_rect1d_all0(EFermi0, p.ems, p.temp, W1, 24e-9, p.nmax, p.mmax)

    fig = plt.figure()
    plt.plot(EFermi,n1,label="8nm, parabollic")
    plt.plot(EFermi,n2,label="16nm, parabollic")
    plt.plot(EFermi,n3,label="24nm, parabollic")
    plt.xlabel('Fermi Energy (eV)')
    plt.ylabel('1DEG Density (cm$^{-1}$)')
    plt.legend(loc='best')
    
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

    EFermi=-0.2
    for i,ene0 in enumerate(ene0_list):
        # val1[i]=func_for_findroot_E0_rect1dNP(ene0,Vds,Vgs,p)
        val1[i]=func_for_findroot_E0_rect1d(ene0,Vds,Vgs,EFermi,p)

    EFermi=0
    for i,ene0 in enumerate(ene0_list):
        # val2[i]=func_for_findroot_E0_rect1dNP(ene0,Vds,Vgs,p)
        val2[i]=func_for_findroot_E0_rect1d(ene0,Vds,Vgs,EFermi,p)

    EFermi=0.2
    for i,ene0 in enumerate(ene0_list):
        # val3[i]=func_for_findroot_E0_rect1dNP(ene0,Vds,Vgs,p)
        val3[i]=func_for_findroot_E0_rect1d(ene0,Vds,Vgs,EFermi,p)
        

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
    
    EFermi=-0.2
    for i, Vgs0 in enumerate(Vgs):
#        val1[i] = E0_rect1dNP_root(Vds,Vgs0,p)
        val1[i] = E0_rect1d_root(Vds,Vgs0,EFermi,p)
#         val1[i] = fetmodel.E0_rect1dNP_root(Vds,Vgs0,p)
    
    EFermi=0
    for i, Vgs0 in enumerate(Vgs):
#        val2[i] = E0_rect1dNP_root(Vds,Vgs0,p)
        val2[i] = E0_rect1d_root(Vds,Vgs0,EFermi,p)
    
    EFermi=0.2
    for i, Vgs0 in enumerate(Vgs):
#        val3[i] = E0_rect1dNP_root(Vds,Vgs0,p)
        val3[i] = E0_rect1d_root(Vds,Vgs0,EFermi,p)

    fig = plt.figure()
    plt.plot(Vgs, val1,label="-0.2")
    plt.plot(Vgs, val2,label="0")
    plt.plot(Vgs, val3,label="0.2")
    plt.xlabel('Gate Voltage (V)')
    plt.ylabel('Top of the barrier height (eV)')
    plt.legend(loc='best')

    ## current
    Vds=0.5
    Vgs = np.linspace(0, 1,endpoint=True)
    val1 = np.empty_like(Vgs)
    val2 = np.empty_like(Vgs)
    val3 = np.empty_like(Vgs)
    
    EFs=-0.2
    for i, Vgs0 in enumerate(Vgs):
        # val1[i] = Ids_ballistic1d_rect1dNP(Vds, Vgs0, p, 0)/(2*(p.W1+p.W2))
        val1[i] = Ids_ballistic1d_rect1d(Vds, Vgs0, p, EFs)/(2*(p.W1+p.W2))
    EFs=0
    for i, Vgs0 in enumerate(Vgs):
        # val2[i] = Ids_ballistic1d_rect1dNP(Vds, Vgs0, p, 0)/(2*(p.W1+p.W2))
        val2[i] = Ids_ballistic1d_rect1d(Vds, Vgs0, p, EFs)/(2*(p.W1+p.W2))
    EFs=0.2
    for i, Vgs0 in enumerate(Vgs):
        # val3[i] = Ids_ballistic1d_rect1dNP(Vds, Vgs0, p, 0)/(2*(p.W1+p.W2))
        val3[i] = Ids_ballistic1d_rect1d(Vds, Vgs0, p, EFs)/(2*(p.W1+p.W2))

    fig, ax = plt.subplots()
    ax.plot(Vgs, val1, label='E$_{Fs}$=-0.2 eV')
    ax.plot(Vgs, val2, label='E$_{Fs}$=0 eV')
    ax.plot(Vgs, val3, label='E$_{Fs}$=0.2 eV')
    ax.set_yscale("log")
    plt.hlines([100e-3], min(Vgs),max(Vgs), "blue", linestyles='dashed')
    plt.xlabel('Gate Voltage (V)')
    plt.ylabel('Drain Current (uA/um)')
    plt.legend(loc='best')
        
    Vds=0.5
    EFs=0
    Vds = np.arange(0, 1, 0.01)
    Ids1 = np.empty_like(Vds)
    Ids2 = np.empty_like(Vds)
    Vgs=0.5
    for i, Vds0 in enumerate(Vds):
        # Ids1[i] = Ids_ballistic1d_rect1dNP(Vds0, Vgs, p, 0)/(2*(p.W1+p.W2))
        Ids1[i] = Ids_ballistic1d_rect1d(Vds0, Vgs, p, EFs)/(2*(p.W1+p.W2))
    Vgs=1
    for i, Vds0 in enumerate(Vds):
        # Ids2[i] = Ids_ballistic1d_rect1dNP(Vds0, Vgs, p, 0)/(2*(p.W1+p.W2))
        # Ids2[i] = fetmodel.Ids_ballistic1d_rect1dNP(Vds0, Vgs, p, 0)
        Ids2[i] = Ids_ballistic1d_rect1d(Vds0, Vgs, p, EFs)/(2*(p.W1+p.W2))

    fig = plt.figure()
    plt.plot(Vds, Ids1, label='Vgs=0.5 V')
    plt.plot(Vds, Ids2, label='Vgs=1 V')
    plt.xlabel('Drain Voltage (V)')
    plt.ylabel('Drain Current (uA/um)')
    plt.legend(loc='best')

    plt.show()
