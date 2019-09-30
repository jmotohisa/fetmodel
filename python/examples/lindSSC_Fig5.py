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
import argparse

# determine EFs to set appropriate Vth
# Ids = 100 nA/um at Vgs=0
def determine_EFs(p, Vds, Ids_cutoff, left=-0.1, right=0.1):
    e0 = optimize.root_scalar(
        func_det_EFs, args=(p, Vds, Ids_cutoff), x0=left, x1=right)
    if(e0.converged==True):
        return e0.root
    else:
        print("EFs convergence error !")
        print(e0)
        return 0


def func_det_EFs(EFs, p, Vds, Ids_cutoff):
    return Ids_cutoff - fetmodel.Ids_ballistic1d_rect1dNP(Vds, 0, p, EFs)/(2*(p.W1+p.W2))


def check_func_det_EFs(p, Vds, Ids_cutoff, left=-0.1, right=0.1):
    ene0_list=np.linspace(left,right,endpoint=True)
    val=np.empty_like(ene0_list)
    for i,ene0 in enumerate(ene0_list):
        val[i]=func_det_EFs(ene0, p, Vds, Ids_cutoff)

    plt.plot(ene0_list,val)
    plt.hlines([0], left, right, "blue", linestyles='dashed')
    plt.xlabel('Source Fermi Energy')
    plt.show()
    return val


def determine_EFs2(p,Vds,Ids_cutoff, left=-0.1, right=0.1):
    e0 = optimize.root_scalar(
        func_det_EFs2, args=(p, Vds, Ids_cutoff), method='brentq',bracket=[-0.1,0.5], x0=left, x1=right)
    if(e0.converged==True):
        return e0.root
    else:
        print("EFs convergence error !")
        print(e0)
        return 0


def func_det_EFs2(EFs, p, Vds, Ids_cutoff):
    retval=Ids_cutoff - ballistic1d.Ids_ballistic1d_rect1dNP(Vds, 0, p, EFs)/(2*(p.W1+p.W2))
    return retval
    

def check_func_det_EFs2(p, Vds, Ids_cutoff, left=-0.1, right=0.1):
    ene0_list=np.linspace(left,right,endpoint=True)
    val=np.empty_like(ene0_list)
    for i,ene0 in enumerate(ene0_list):
        val[i]=func_det_EFs2(ene0, p, Vds, Ids_cutoff)

    plt.plot(ene0_list,val)
    plt.hlines([0], left, right, "blue", linestyles='dashed')
    plt.xlabel('Source Fermi Energy')
    plt.show()
    return val

def get_args():
    # 準備
    parser = argparse.ArgumentParser(
        description='Rerpoduce Figs. 5(a) and (b) in E.Lind')
    
    # 標準入力以外の場合
    parser.add_argument('-n', '--nmax',
                        nargs='?',
                        type=int,
                        help='nmax',
                        default=5)
    parser.add_argument('-m', '--mmax',
                        nargs='?',
                        type=int,
                        help='mmax',
                        default=5)
    # parser.add_argument("--alert", help="optional", action="store_true")
    
    # 結果を受ける
    args = parser.parse_args()
    
    return(args)


if __name__=='__main__':
    args = get_args()
    nmax = args.nmax
    mmax = args.mmax

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
    p=fetmodel.parameters_ballistic(alpha=alpha,
                                    Ceff=Cox*Cc/(Cox+Cc),
                                    ems=ems,
                                    W1=W1,
                                    W2=W2,
                                    nmax=nmax,
                                    mmax=mmax)
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
    EFs = determine_EFs(p, Vds, Ids_cutoff,left=0., right=0.1)
    print("Source Fermi Energy for appropriate Vth (eV):", EFs)
    Ids0 = fetmodel.Ids_ballistic1d_rect1dNP(Vds, Vgs, p, EFs)/(2*(p.W1+p.W2))*1e3
    print("Actual Ids at Vgs=0 (nA/um): ", Ids0)
    print("")

    check_func_det_EFs2(p, Vds, Ids_cutoff)
    print("Determine EFs using ballistic1d")
    EFs2 = determine_EFs2(p, Vds, Ids_cutoff,left=0., right=0.1) # ,left=0.2,right=0.3
    print("Source Fermi Energy for appropriate Vth (eV):", EFs2)
    Ids0 = ballistic1d.Ids_ballistic1d_rect1dNP(Vds, Vgs, p, EFs2)/(2*(p.W1+p.W2))*1e3
    print("Actual Ids at Vgs=0 (nA/um): ", Ids0)
    print("")

    dVgs = 0.005
    Vgs = np.arange(0, 0.6, dVgs)
    Ids1 = np.empty_like(Vgs)
    Ids2 = np.empty_like(Vgs)
    Ids3 = np.empty_like(Vgs)

    p.nmax=2
    p.mmax=1
    for i, Vgs0 in enumerate(Vgs):
        # Ids1[i] = fetmodel.Ids_ballistic1d_rect1dNP(
        #     Vds, Vgs0, p, EFs)/(2*(p.W1+p.W2))*1e-3
        Ids1[i] = ballistic1d.Ids_ballistic1d_rect1dNP(
            Vds, Vgs0, p, EFs2)/(2*(p.W1+p.W2))*1e-3
    lstr1='(n,m)=('+str(p.nmax)+','+str(p.mmax)+')'

    p.nmax=2
    p.mmax=2
    for i, Vgs0 in enumerate(Vgs):
        # Ids1[i] = fetmodel.Ids_ballistic1d_rect1dNP(
        #     Vds, Vgs0, p, EFs)/(2*(p.W1+p.W2))*1e-3
        Ids2[i] = ballistic1d.Ids_ballistic1d_rect1dNP(
            Vds, Vgs0, p, EFs2)/(2*(p.W1+p.W2))*1e-3
    lstr2='(n,m)=('+str(p.nmax)+','+str(p.mmax)+')'

    p.nmax=nmax
    p.mmax=mmax
    for i, Vgs0 in enumerate(Vgs):
        # Ids1[i] = fetmodel.Ids_ballistic1d_rect1dNP(
        #     Vds, Vgs0, p, EFs)/(2*(p.W1+p.W2))*1e-3
        Ids3[i] = ballistic1d.Ids_ballistic1d_rect1dNP(
            Vds, Vgs0, p, EFs2)/(2*(p.W1+p.W2))*1e-3
    lstr3='(n,m)=('+str(p.nmax)+','+str(p.mmax)+')'

    gm1 = np.gradient(Ids1, dVgs)
    gm2 = np.gradient(Ids2, dVgs)
    gm3 = np.gradient(Ids3, dVgs)

    plt.figure(figsize=(4,6))
    plt.plot(Vgs, Ids1, label=lstr1)
    plt.plot(Vgs, Ids2, label=lstr2)
    plt.plot(Vgs, Ids3, label=lstr3)
    plt.xlabel('Vgs - Vth (V)',fontsize=18)
    plt.ylabel('Ids (mA/um)',fontsize=18)
    plt.legend(loc='best',fontsize=18)
    plt.title('Vds='+str(Vds)+' V',fontsize=18)
    plt.xlim([0,0.6])
    plt.ylim([0,1.5])
    plt.tick_params(labelsize=18)
    plt.tight_layout()

    fig, ax = plt.subplots(figsize=(4,6))
        # "fig1=plt.figure(figsize=(6,4))\n",
        #     "ax1=fig1.add_subplot(111)\n",
    ax.plot(Vgs, Ids1, label=lstr1)
    ax.plot(Vgs, Ids2, label=lstr2)
    ax.plot(Vgs, Ids3, label=lstr3)
    ax.set_yscale("log")
    plt.xlabel('Vgs - Vth (V)', fontsize=18)
    plt.ylabel('Ids (mA/um)', fontsize=18)
    plt.legend(loc='best', fontsize=18)
    plt.hlines([Ids_cutoff*1e-3], min(Vgs),max(Vgs), "blue", linestyles='dashed')
    plt.title('Vds='+str(Vds)+' V', fontsize=18)
    plt.tick_params(labelsize=18)
    plt.tight_layout()

    ## transconductance
    plt.figure(figsize=(4,6))
    plt.plot(Vgs, gm1, label=lstr1)
    plt.plot(Vgs, gm2, label=lstr2)
    plt.plot(Vgs, gm3, label=lstr3)
    plt.xlabel('Vgs - Vth (V)', fontsize=18)
    plt.ylabel('Transconductance (mS/um)', fontsize=18)
    plt.legend(loc='best', fontsize=18)
    plt.title('Vds='+str(Vds)+' V', fontsize=18)
    plt.tick_params(labelsize=18)
    plt.xlim([0,0.6])
    plt.ylim([0,5])
    plt.tight_layout()

    plt.show()

    for i, Vgs0 in enumerate(Vgs):
        print(Vgs0,Ids1[i],Ids2[i],Ids3[i],gm1[i],gm2[i],gm3[i])
