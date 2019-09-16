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
from lindSSC_Fig5 import *
import pdb

def func_ems(Eg,q20=20):
    return 1/(1+q20/Eg)

def set_p(p,Eg,W2,nmax,mmax):
    p.ems = func_ems(Eg)
    p.alpha = fetmodel.alpha_NP(Eg, p.ems)
    p.W2 = W2
    Cox = fetmodel.Cox_rect(epsOX, tOX, p.W1, W2)
    Cc = fetmodel.Cc_rect(epsS, p.W1, W2)
    p.Ceff = Cox*Cc/(Cox+Cc)
    p.nmax = nmax
    p.mmax = mmax

def IdsVgs(p,Vds,Vgs_list,EFs):
    Ids=np.empty_like(Vgs_list)
    for i,Vgs0 in enumerate(Vgs_list):
        Ids[i]=ballistic1d.Ids_ballistic1d_rect1dNP(Vds, Vgs0, p, EFs)/(2*(p.W1+p.W2))
    return Ids
    

if __name__ == '__main__':
    EFermi=0
    Eg = 0.36
    epsOX = 20
    epsS = 15.15
    tOX = 3e-9
    temperature = 300
    ems = 0.023
    W1 = 5e-9
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
    Ids_cutoff=100e-9*1e6
    
    # Figure 6:
    Vds=0.5
    Vgs=0.5
    Eg_list=np.linspace(0.3,3.4,endpoint=True)
    Ion1=np.empty_like(Eg_list)
    Ion2=np.empty_like(Eg_list)
    Ion3=np.empty_like(Eg_list)
    nmax=3
    mmax=3
    for i, Eg in enumerate(Eg_list):
        W2=8e-9
        set_p(p,Eg,W2,nmax,mmax)
        # EFs = fetmodel.determine_EFs(p, Vds, Ids_cutoff)
        # Ion1[i] = fetmodel.Ids_ballistic1d_rect1dNP(Vds, Vgs, p, EFs)/(2*(p.W1+p.W2))*1e-3
        EFs2 = determine_EFs2(p, Vds, Ids_cutoff)
        # print('Ioff(nA/mm)',ballistic1d.Ids_ballistic1d_rect1dNP(Vds,0 , p, EFs2)/(2*(p.W1+p.W2))*1e3)
        Ion1[i] = ballistic1d.Ids_ballistic1d_rect1dNP(Vds, Vgs, p, EFs2)/(2*(p.W1+p.W2))*1e-3

        W2=16e-9
        set_p(p,Eg,W2,nmax,mmax)
        # EFs = fetmodel.determine_EFs(p, Vds, Ids_cutoff)
        # Ion2[i] = fetmodel.Ids_ballistic1d_rect1dNP(Vds, Vgs, p, EFs)/(2*(p.W1+p.W2))*1e-3
        EFs2 = determine_EFs2(p, Vds, Ids_cutoff)
        # print('Ioff(nA/um)=',ballistic1d.Ids_ballistic1d_rect1dNP(Vds,0 , p, EFs2)/(2*(p.W1+p.W2))*1e3)
        Ion2[i] = ballistic1d.Ids_ballistic1d_rect1dNP(Vds, Vgs, p, EFs2)/(2*(p.W1+p.W2))*1e-3

        W2=24e-9
        set_p(p,Eg,W2,nmax,mmax)
        # EFs = fetmodel.determine_EFs(p, Vds, Ids_cutoff)
        # Ion3[i] = fetmodel.Ids_ballistic1d_rect1dNP(Vds, Vgs, p, EFs)/(2*(p.W1+p.W2))*1e-3
        EFs2 = determine_EFs2(p, Vds, Ids_cutoff)
        # print('Ioff(nA/um)=',ballistic1d.Ids_ballistic1d_rect1dNP(Vds,0 , p, EFs2)/(2*(p.W1+p.W2))*1e3)
        Ion3[i] = ballistic1d.Ids_ballistic1d_rect1dNP(Vds, Vgs, p, EFs2)/(2*(p.W1+p.W2))*1e-3
        
    plt.plot(Eg_list, Ion1, label='8nm')
    plt.plot(Eg_list, Ion2, label='16nm')
    plt.plot(Eg_list, Ion3, label='24nm')
    plt.xlabel('$E_G$ (eV)',fontsize=18)
    plt.ylabel('Ion (mA/um)',fontsize=18)
    plt.legend(loc='best',fontsize=18)
    plt.tick_params(labelsize=18)
    plt.ylim([0.4,1.4])
    plt.title('(n,m)=('+str(p.nmax)+','+str(p.mmax)+')',fontsize=18)
    plt.tight_layout()
    plt.show()

    ####
    W2=24e-9
    Vds=0.5
    Vgs=np.linspace(-0.2,2,endpoint=True)

    Eg=3.4
    set_p(p,Eg,W2,3,3)
    p.output()
    EFs2 = determine_EFs2(p, Vds, Ids_cutoff)
    Ids1=IdsVgs(p,Vds,Vgs,EFs2)
    set_p(p,Eg,W2,5,5)
    EFs2 = determine_EFs2(p, Vds, Ids_cutoff)
    Ids2=IdsVgs(p,Vds,Vgs,EFs2)

    fig, ax = plt.subplots()
    ax.plot(Vgs, Ids1, label='3,3')
    ax.plot(Vgs, Ids2, label='4,4')
    ax.set_yscale("log")
    plt.xlabel('Gate Voltage Vgs-Vth',fontsize=18)
    plt.ylabel('Ids (uA/um)',fontsize=18)
    plt.legend(loc='best',fontsize=18)
    plt.tick_params(labelsize=18)
    plt.title('Eg='+str(Eg)+' eV',fontsize=18)
    plt.tight_layout()

    plt.show()
