#!/usr/bin/env python
# -*- coding: utf-8 -*-

import fetmodel
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import scipy.constants as const

if __name__ == '__main__':
    # Eg = 3.4
    # ems = 0.2
    # epsOX = 8.5
    # epsS = 8.9
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

    print("Test of density1d: parabolic and nonparabolic band")

    ### energy levels in parabolic and nonparabollic bands
    W2_list=np.linspace(5e-9,30e-9,endpoint=True)
    Enmp11=np.empty_like(W2_list)
    Enmp12=np.empty_like(W2_list)
    Enmp21=np.empty_like(W2_list)
    Enmp22=np.empty_like(W2_list)
    Enm11=np.empty_like(W2_list)
    Enm12=np.empty_like(W2_list)
    Enm21=np.empty_like(W2_list)
    Enm22=np.empty_like(W2_list)
    for i,W2 in enumerate(W2_list):
        Enmp11[i] = fetmodel.Ep_nm_rect1d(p.ems, p.W1, W2, 1,1)
        gamma_nm = fetmodel.gamma_nm_NP(Enmp11[i], p.alpha)
        Enm11[i] = fetmodel.E_nm_NP(p.alpha, gamma_nm)
        Enmp12[i] = fetmodel.Ep_nm_rect1d(p.ems, p.W1, W2, 1,2)
        gamma_nm = fetmodel.gamma_nm_NP(Enmp12[i], p.alpha)
        Enm12[i] = fetmodel.E_nm_NP(p.alpha, gamma_nm)
        Enmp21[i] = fetmodel.Ep_nm_rect1d(p.ems, p.W1, W2, 2,1)
        gamma_nm = fetmodel.gamma_nm_NP(Enmp21[i], p.alpha)
        Enm21[i] = fetmodel.E_nm_NP(p.alpha, gamma_nm)
        Enmp22[i] = fetmodel.Ep_nm_rect1d(p.ems, p.W1, W2, 2,2)
        gamma_nm = fetmodel.gamma_nm_NP(Enmp22[i], p.alpha)
        Enm22[i] = fetmodel.E_nm_NP(p.alpha, gamma_nm)

    fig = plt.figure()
    plt.plot(W2_list,Enmp11,label='(n,m)=(1,1), parabolic',color='black')
    plt.plot(W2_list,Enm11,label='(n,m)=(1,1), nonparabolic',color='black',linestyle='dashed')
    plt.plot(W2_list,Enmp12,label='(n,m)=(1,2), parabolic',color='red')
    plt.plot(W2_list,Enm12,label='(n,m)=(1,2), nonparabolic',color='red',linestyle='dashed')
    plt.plot(W2_list,Enmp21,label='(n,m)=(2,1), parabolic',color='magenta')
    plt.plot(W2_list,Enm21,label='(n,m)=(2,1), nonparabolic',color='magenta',linestyle='dashed')
    plt.plot(W2_list,Enmp22,label='(n,m)=(2,2), parabolic',color='blue')
    plt.plot(W2_list,Enm22,label='(n,m)=(2,2), nonparabolic',color='blue',linestyle='dashed')
    plt.xlabel('width W2 (m)')
    plt.ylabel('Energy (eV)')
    plt.legend(loc='best')
    plt.title('Energy Levels in a NW: W1='+str(p.W1)+' m')
    plt.show()

    ### 1DEG density
    EFermi=np.linspace(-0.1,0.5,endpoint=True)
    n1p=np.empty_like(EFermi)
    n2p=np.empty_like(EFermi)
    n3p=np.empty_like(EFermi)
    n1=np.empty_like(EFermi)
    n2=np.empty_like(EFermi)
    n3=np.empty_like(EFermi)

    # pdb.set_trace()
    ## parabolic band
    for i,EFermi0 in enumerate(EFermi):
        n1p[i]=fetmodel.density1d_rect1d_all0(EFermi0, p.ems, p.temp, W1, 8e-9, p.nmax, p.mmax)
        n2p[i]=fetmodel.density1d_rect1d_all0(EFermi0, p.ems, p.temp, W1, 16e-9, p.nmax, p.mmax)
        n3p[i]=fetmodel.density1d_rect1d_all0(EFermi0, p.ems, p.temp, W1, 24e-9, p.nmax, p.mmax)
        
    ## nonparabolic band
    for i,EFermi0 in enumerate(EFermi):
        n1[i]=fetmodel.density1d_rect1dNP_all0(
            EFermi0, p.alpha, p.ems, p.temp, W1, 8e-9, p.nmax, p.mmax)
        n2[i]=fetmodel.density1d_rect1dNP_all0(
            EFermi0, p.alpha, p.ems, p.temp, W1, 16e-9, p.nmax, p.mmax)
        n3[i]=fetmodel.density1d_rect1dNP_all0(
            EFermi0, p.alpha, p.ems, p.temp, W1, 24e-9, p.nmax, p.mmax)

    fig = plt.figure()
    plt.plot(EFermi,n1p,label="8nm, parabolic", linestyle='dashed',color='red')
    plt.plot(EFermi,n2p,label="16nm, parabolic", linestyle='dashed',color='black')
    plt.plot(EFermi,n3p,label="24nm, parabolic", linestyle='dashed',color='blue')
    plt.plot(EFermi,n1,label="8nm, nonparabolic",color='red')
    plt.plot(EFermi,n2,label="16nm, nonparabolic",color='black')
    plt.plot(EFermi,n3,label="24nm, nonparabolic",color='blue')
    plt.xlabel('Fermi Energy (eV)')
    plt.ylabel('1DEG Density (m$^{-1}$)')
    plt.xlim([-0.1,0.2])
    # plt.ylim([0,1.5e9])
    plt.legend(loc='best')

    fig, ax = plt.subplots()
    ax.plot(EFermi,n1p,label="8nm, parabolic", linestyle='dashed',color='red')
    ax.plot(EFermi,n2p,label="16nm, parabolic", linestyle='dashed',color='black')
    ax.plot(EFermi,n3p,label="24nm, parabolic", linestyle='dashed',color='blue')
    ax.plot(EFermi,n1,label="8nm, nonparabolic",color='red')
    ax.plot(EFermi,n2,label="16nm, nonparabolic",color='black')
    ax.plot(EFermi,n3,label="24nm, nonparabolic",color='blue')
    plt.xlabel('Fermi Energy (eV)')
    plt.ylabel('1DEG Density (m$^{-1}$)')
    plt.legend(loc='best')
    ax.set_yscale("log")

    plt.show()

