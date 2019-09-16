#!/usr/bin/env python
# -*- coding: utf-8 -*-

# From : Lind, E. (2016). High frequency III{\textendash}V nanowire MOSFETs.
# Semiconductor Science and Technology, 31(9), 93005–93014.
# https://doi.org/10.1088/0268-1242/31/9/093005

# Fig. 6(a)

import fetmodel
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import scipy.constants as const
import ballistic1d
import lindSSC_Fig5
import pdb

def func_ems(Eg,q20=20):
    return 1/(1+20/Eg)

def set_p(p,Eg,W2,nmax,mmax):
    p.ems = func_ems(Eg)
    p.alpha = fetmodel.alpha_NP(Eg, p.ems)
    # p.W1 = W1
    p.W2 = W2
    # p.EFermi = 0
    # p.alpha_D = alpha_D
    # p.alpha_G = alpha_G
    Cox = fetmodel.Cox_rect(epsOX, tOX, p.W1, W2)
    Cc = fetmodel.Cc_rect(epsS, p.W1, W2)
    p.Ceff = Cox*Cc/(Cox+Cc)
    # p.temp = temperature
    p.nmax = nmax
    p.mmax = mmax


"""
determine EFs to set appropriate Vth
Ids = 100 nA/um at Vgs=0
"""
def determine_EFs(p, Vds, Ids):
    e0 = optimize.root_scalar(
        func_det_EFs, args=(p, Vds, Ids), x0=-0.1, x1=0.1)
    if(e0.converged==True):
        return e0.root
    else:
        print("EFs convergence error !")
        print(e0)
        return 0


def func_det_EFs(EFs, p, Vds, Ids):
    return Ids - fetmodel.Ids_ballistic1d_rect1dNP(Vds, 0, p, EFs)/(2*(p.W1+p.W2))


if __name__ == '__main__':
    Eg = 0.36
    epsOX = 20
    epsS = 15.15
    tOX = 3e-9
    temperature = 300
    ems = 0.023
    W1 = 5e-9
    W2 = 8e-9
    alpha = fetmodel.alpha_NP(Eg, ems)
    alpha_D = 0
    alpha_G = 1
    p = fetmodel.param_ballistic_new()

    Vds = 0.5
    Vgs = 0.5
    nmax=2
    mmax=2
    set_p(p,Eg,8e-9,nmax,mmax)
    print(p.ems,p.alpha)

    Ids_cutoff=100e-9*1e6
    # EFs_list = np.arange(-0.1,0.4,0.005)
    # funcval=np.empty_like(EFs_list)
    # for i,EFs in enumerate(EFs_list):
    #     funcval[i] = ballistic_lindSSC1.func_det_EFs2(EFs, p, Vds, Ids_cutoff)

    # plt.plot(EFs_list,funcval)
    # plt.show()
    
    left=0.25
    right=0.35
    EFs = ballistic_lindSSC1.determine_EFs2(p, Vds, Ids_cutoff,left,right)
    print("Source Fermi Energy for appropriate Vth (eV):", EFs)
    Ids0 = fetmodel.Ids_ballistic1d_rect1dNP(Vds, 0, p, EFs)/(2*(p.W1+p.W2))*1e3
    print("Actual Ids at Vgs=0 (nA/um): ", Ids0)

    # Figure 6:
    Eg_list=np.arange(0.3,3.2,0.05)
    Ion1=np.empty_like(Eg_list)
    Ion2=np.empty_like(Eg_list)
    Ion3=np.empty_like(Eg_list)
    nmax=3
    mmax=3
    for i, Eg in enumerate(Eg_list):
        set_p(p,Eg,8e-9,nmax,mmax)
        # pdb.set_trace()
        EFs = ballistic_lindSSC1.determine_EFs2(p, Vds, Ids_cutoff,left,right)
        # EFs = determine_EFs(p, Vds, Ids_cutoff)
        # Ion1[i] = fetmodel.Ids_ballistic1d_rect1dNP(Vds, Vgs, p, EFs)/(2*(p.W1+p.W2))*1e-3
        Ion1[i] = ballistic1d.func_current1D(Vgs, Vds, p, EFs)/(2*(p.W1+p.W2))*1e-3

        set_p(p,Eg,16e-9,nmax,mmax)
        EFs = ballistic_lindSSC1.determine_EFs2(p, Vds, Ids_cutoff,left,right)
        # EFs = determine_EFs(p, Vds, Ids_cutoff)
        # Ion1[i] = fetmodel.Ids_ballistic1d_rect1dNP(Vds, Vgs, p, EFs)/(2*(p.W1+p.W2))*1e-3
        Ion2[i] = ballistic1d.func_current1D(Vgs, Vds, p, EFs)/(2*(p.W1+p.W2))*1e-3

        set_p(p,Eg,24e-9,nmax,mmax)
        EFs = ballistic_lindSSC1.determine_EFs2(p, Vds, Ids_cutoff,left,right)
        # EFs = determine_EFs(p, Vds, Ids_cutoff)
        # Ion1[i] = fetmodel.Ids_ballistic1d_rect1dNP(Vds, Vgs, p, EFs)/(2*(p.W1+p.W2))*1e-3
        Ion3[i] = ballistic1d.func_current1D(Vgs, Vds, p, EFs)/(2*(p.W1+p.W2))*1e-3
        
        
# dVgs = 0.005
# Vgs = np.arange(0, 0.6, dVgs
# Ids1 = np.empty_like(Vgs)
# for i, Vgs0 in enumerate(Vgs):
#     Ids1[i] = fetmodel.Ids_ballistic1d_rect1dNP(
#         Vds, Vgs0, p, EFs)/(2*(p.W1+p.W2))*1e-3

# gm1 = np.gradient(Ids1, dVgs)

    plt.plot(Eg_list, Ion1, label='8nm')
    plt.plot(Eg_list, Ion2, label='16nm')
    plt.plot(Eg_list, Ion3, label='24nm')
    plt.xlabel('$E_G$ (eV)')
    plt.ylabel('Ion (mA/um)')
    plt.legend(loc='best')
    plt.show()

    # plt.plot(Vgs, gm1, label='gm')
    # plt.xlabel('Vgs - Vth (V)')
    # plt.ylabel('Transconductance (mS/um)')
    # plt.show()
