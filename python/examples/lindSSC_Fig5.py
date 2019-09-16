#!/usr/bin/env python
# -*- coding: utf-8 -*-

# From : Lind, E. (2016). High frequency III{\textendash}V nanowire MOSFETs.
# Semiconductor Science and Technology, 31(9), 93005–93014.
# https://doi.org/10.1088/0268-1242/31/9/093005
# Figure 5:

import fetmodel
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import ballistic1d

# determine EFs to set appropriate Vth
# Ids = 100 nA/um at Vgs=0
def determine_EFs(p, Vds, Ids_cutoff, left=-0.1, right=0.1):
    e0 = optimize.root_scalar(
        func_det_EFs, args=(p, Vds, Ids_cutoff), x0=left, x1=right)
    return e0.root


def func_det_EFs(EFs, p, Vds, Ids_cutoff):
    return Ids_cutoff - fetmodel.Ids_ballistic1d_rect1dNP(Vds, 0, p, EFs)/(2*(p.W1+p.W2))


def check_func_det_EFs(p, Vds, Ids_cutoff, left=0.2, right=0.3):
    ene0_list=np.linspace(left,right,endpoint=True)
    val=np.empty_like(ene0_list)
    for i,ene0 in enumerate(ene0_list):
        val[i]=func_det_EFs(ene0, p, Vds, Ids_cutoff)

    plt.plot(ene0_list,val)
    plt.hlines([0], left, right, "blue", linestyles='dashed')
    plt.xlabel('Source Fermi Energy')
    plt.show()
    return val


def determine_EFs2(p,Vds,Ids_cutoff, left=0.2, right=0.3):
    e0 = optimize.root_scalar(
        func_det_EFs2, args=(p, Vds, Ids_cutoff), x0=left, x1=right)
    if(e0.converged==True):
        return e0.root
    else:
        print("EFs convergence error !")
        print(e0)
        return 0


def func_det_EFs2(EFs, p, Vds, Ids_cutoff):
    retval=Ids_cutoff - ballistic1d.Ids_ballistic1d_rect1dNP(Vds, 0, p, EFs)/(2*(p.W1+p.W2))
    return retval
    

def check_func_det_EFs2(p, Vds, Ids_cutoff, left=0.2, right=0.3):
    ene0_list=np.linspace(left,right,endpoint=True)
    val=np.empty_like(ene0_list)
    for i,ene0 in enumerate(ene0_list):
        val[i]=func_det_EFs2(ene0, p, Vds, Ids_cutoff)

    plt.plot(ene0_list,val)
    plt.hlines([0], left, right, "blue", linestyles='dashed')
    plt.xlabel('Source Fermi Energy')
    plt.show()
    return val

if __name__=='__main__':
    EFermi=-0.1
    Eg = 0.36
    epsOX = 20
    epsS = 15.15
    tOX = 3e-9
    temperature = 300
    ems = 0.023
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
    print("")
    
    # p = fetmodel.param_ballistic_new()
    # p.ems = ems
    # p.alpha = alpha
    # p.W1 = W1
    # p.W2 = W2
    # p.EFermi = 0
    # p.alpha_D = alpha_D
    # p.alpha_G = alpha_G
    # p.Ceff = Cox*Cc/(Cox+Cc)
    # p.temp = temperature
    # p.nmax = 3
    # p.mmax = 3
    
    Vds = 0.5
    Vgs = 0
    
    # EFs = np.arange(-0.1, 0.1, 0.005)
    # Ids1 = np.empty_like(EFs)
    # for i, EFs0 in enumerate(EFs):
    #     Ids1[i] = fetmodel.Ids_ballistic1d_rect1dNP(
    #         Vds, Vgs, p, EFs0)/(2*(p.W1+p.W2))*1e3
    
    # plt.plot(EFs, Ids1, label='Ids')
    # plt.xlabel('Source Fermi Energy (eV)')
    # plt.ylabel('Drain Current (nA/um)')
    # plt.show()
    
    # print(fetmodel.Ids_ballistic1d_rect1dNP(Vds, Vgs, p, 0.075)/(2*(p.W1+p.W2))*1e3)
    # ballistic1d.check_func_for_E0_rect1dNP_fetmodel(-0.2,0.2,Vds,Vgs,p)

    Ids_cutoff=100e-9*1e6

    check_func_det_EFs(p, Vds, Ids_cutoff)
    print("Determine EFs using fetmodel")
    EFs = determine_EFs(p, Vds, Ids_cutoff)
    print("Source Fermi Energy for appropriate Vth (eV):", EFs)
    Ids0 = fetmodel.Ids_ballistic1d_rect1dNP(Vds, Vgs, p, EFs)/(2*(p.W1+p.W2))*1e3
    print("Actual Ids at Vgs=0 (nA/um): ", Ids0)
    print("")

    check_func_det_EFs2(p, Vds, Ids_cutoff)
    print("Determine EFs using ballistic1d")
    EFs2 = determine_EFs2(p, Vds, Ids_cutoff) # ,left=0.2,right=0.3
    print("Source Fermi Energy for appropriate Vth (eV):", EFs2)
    Ids0 = ballistic1d.Ids_ballistic1d_rect1dNP(Vds, Vgs, p, EFs2)/(2*(p.W1+p.W2))*1e3
    print("Actual Ids at Vgs=0 (nA/um): ", Ids0)
    print("")

    dVgs = 0.005
    Vgs = np.arange(0, 0.6, dVgs)
    Ids1 = np.empty_like(Vgs)
    for i, Vgs0 in enumerate(Vgs):
        # Ids1[i] = fetmodel.Ids_ballistic1d_rect1dNP(
        #     Vds, Vgs0, p, EFs)/(2*(p.W1+p.W2))*1e-3
        Ids1[i] = ballistic1d.Ids_ballistic1d_rect1dNP(
            Vds, Vgs0, p, EFs2)/(2*(p.W1+p.W2))*1e-3
        
    gm1 = np.gradient(Ids1, dVgs)

    plt.plot(Vgs, Ids1, label='Ids1')
    plt.xlabel('Vgs - Vth (V)')
    plt.ylabel('Ids (mA/um)')
    plt.show()
    
    plt.plot(Vgs, gm1, label='gm')
    plt.xlabel('Vgs - Vth (V)')
    plt.ylabel('Transconductance (mS/um)')
    plt.show()
