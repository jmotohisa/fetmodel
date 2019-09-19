#!/usr/bin/env python
# coding: utf-8

import fetmodel
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import scipy.constants as const
# from scipy import integrate
# import pdb

####
### Parabolic band
###
def density1d_circ1d_all0(EFermi, ems, temp, radius, nmax, mmax):
    """ electron denstiy in circular nanowire (parabolic band)
    """
    n0 = 0
    nlist = np.arange(1, nmax+1, dtype=np.int64)
    mlist = np.arange(1, mmax+1, dtype=np.int64)
    for n in nlist:
        for m in mlist:
            Enmp = fetmodel.Ep_nm_radial1d(ems, radius, int(n), int(m))
            # print(Enmp)
            nn = fetmodel.density1d0(EFermi, Enmp, ems, temp)
            n0 += nn

    return(n0)

def density1d_circ1d_all(p):
    return density1d_circ1d_all0(p.EFermi, p.ems, p.temp, p.W1/2, p.nmax, p.mmax)
    

def func_for_findroot_E0_circ1d(ene0, Vds, Vgs, p):
    """ function to find find energy of the top of the barrier in circular NW
    """
    n1d_S = density1d_circ1d_all0(p.EFermi - ene0, p.ems, p.temp, p.W1/2, p.nmax,p.mmax)
    n1d_D = density1d_circ1d_all0(
        p.EFermi - ene0 - Vds, p.ems, p.temp, p.W1/2, p.nmax,p.mmax)
    q0 = const.elementary_charge * (n1d_S + n1d_D) / (2 * p.Ceff)
    return ene0 + (p.alpha_D * Vds + p.alpha_G * Vgs - q0)


def check_func_for_E0_circ1d(Vds,Vgs,p,left=-0.2,right=0):
    """
    Plot function (Python implementation) to find E0
    Parabolic band
    """
    left0 = -(p.alpha_D*Vds+p.alpha_G*Vgs) - 0.2
    left2=min([left,left0])
    right0=-(p.alpha_D*Vds+p.alpha_G*Vgs)
    right2=max([right,right0])
    ene0_list=np.linspace(left2,right2,endpoint=True)
    val=np.empty_like(ene0_list)
    for i,ene0 in enumerate(ene0_list):
        val[i]=func_for_findroot_E0_circ1d(ene0, Vds, Vgs, p)

    plt.plot(ene0_list,val)
    plt.show()
    return val


def E0_circ1d_root(Vds,Vgs,p,left=-0.2,right=0):
    """ get energy of the top of the barrier
    """
    left0 = -(p.alpha_D*Vds+p.alpha_G*Vgs) - 0.2
    left2=min([left,left0])
    right0=-(p.alpha_D*Vds+p.alpha_G*Vgs)
    right2=max([right,right0])

    e0 = optimize.root_scalar(
        func_for_findroot_E0_circ1d, args=(Vds, Vgs, p), x0=-0.0, x1=0.15)
    if(e0.converged==True):
        return e0.root
    else:
        print("E0 convergence error !")
        print(e0)
        return 0


def Ids_ballistic1d_circ1d(Vds, Vgs, p, EFs,left=-0.2,right=0):
    """ current of ballistic FET with circular 1d NW
    """
    e0=E0_circ1d_root(Vds,Vgs,p,left,right)
    # e00 = optimize.root_scalar(
    #     func_e0_find_circ1d, args=(Vgs, Vds, p), x0=-0.1, x1=1)
    # e0=e00.root
    nlist = np.arange(1, p.nmax+1, dtype=np.int64)
    mlist = np.arange(1, p.mmax+1, dtype=np.int64)
    cur = 0
    for m in mlist:
        for n in nlist:
            Enp = fetmodel.Ep_nm_radial1d(p.ems, p.W1/2, int(n), int(m))
            cur1 = func_FD0(EFs-Enp-e0, p.temp)
            cur2 = func_FD0(EFs-Enp-e0-Vds, p.temp)
            cur += cur1-cur2

    return (cur*2*const.elementary_charge/const.h*p.temp*const.Boltzmann)


####
### Nonparabolic band
###
def density1d_circ1dNP_all0(EFermi, alpha, ems, temp, radius, nmax, mmax):
    """ electron denstiy in circular nanowire (parabolic band)
    """
    n0 = 0
    nlist = np.arange(1, nmax+1, dtype=np.int64)
    mlist = np.arange(1, mmax+1, dtype=np.int64)
    for n in nlist:
        for m in mlist:
            Enmp = fetmodel.Ep_nm_radial1d(ems, radius, int(n), int(m))
            gamma_nm = fetmodel.gamma_nm_NP(Enmp, alpha)
            alpha_nm = fetmodel.alpha_nm_NP(alpha,gamma_nm)
            ems_nm = fetmodel.ems_nm_NP(ems,gamma_nm)
            Enm = fetmodel.E_nm_NP(alpha, gamma_nm)
            # print(Enm)
            nn = fetmodel.density1d_NP0(EFermi, Enm, alpha_nm, ems_nm, temp)
            n0 += nn

    return(n0)

def density1d_circ1dNP_all(p):
    return density1d_circ1d_all0(p.EFermi, p.alpha, p.ems, p.temp, p.W1/2, p.nmax, p.mmax)
    

def func_for_findroot_E0_circ1dNP(ene0, Vds, Vgs, p):
    """ function to find find energy of the top of the barrier in circular NW
    """
    n1d_S = density1d_circ1dNP_all0(p.EFermi - ene0, p.alpha, p.ems, p.temp, p.W1, p.nmax,p.mmax)
    n1d_D = density1d_circ1dNP_all0(
        p.EFermi - ene0 - Vds, p.alpha, p.ems, p.temp, p.W1, p.nmax,p.mmax)
    q0 = const.elementary_charge * (n1d_S + n1d_D) / (2 * p.Ceff)
    return ene0 + (p.alpha_D * Vds + p.alpha_G * Vgs - q0)


def check_func_for_E0_circ1dNP(Vds,Vgs,p,left=-0.2,right=0):
    """
    Plot function (Python implementation) to find E0
    Parabolic band
    """
    left0 = -(p.alpha_D*Vds+p.alpha_G*Vgs) - 0.2
    left2=min([left,left0])
    right0=-(p.alpha_D*Vds+p.alpha_G*Vgs)
    right2=max([right,right0])
    ene0_list=np.linspace(left2,right2,endpoint=True)
    val=np.empty_like(ene0_list)
    for i,ene0 in enumerate(ene0_list):
        val[i]=func_for_findroot_E0_circ1dNP(ene0, Vds, Vgs, p)

    plt.plot(ene0_list,val)
    plt.show()
    return val


def E0_circ1dNP_root(Vds,Vgs,p,left=-0.2,right=0):
    """ get energy of the top of the barrier
    """
    left0 = -(p.alpha_D*Vds+p.alpha_G*Vgs) - 0.2
    left2=min([left,left0])
    right0=-(p.alpha_D*Vds+p.alpha_G*Vgs)
    right2=max([right,right0])

    e0 = optimize.root_scalar(
        func_for_findroot_E0_circ1dNP, args=(Vgs, Vds, p), x0=-0.0, x1=0.15)
    if(e0.converged==True):
        return e0.root
    else:
        print("E0 convergence error !")
        print(e0)
        return 0


def Ids_ballistic1d_circ1dNP(Vds, Vgs, p, EFs,left=-0.2,right=0):
    """ current of ballistic FET with circular 1d NW
    """
    e0=E0_circ1dNP_root(Vds,Vgs,p,left,right)
    # e00 = optimize.root_scalar(
    #     func_e0_find_circ1d, args=(Vgs, Vds, p), x0=-0.1, x1=1)
    # e0=e00.root
    nlist = np.arange(1, p.nmax+1, dtype=np.int64)
    mlist = np.arange(1, p.mmax+1, dtype=np.int64)
    cur = 0
    for m in mlist:
        for n in nlist:
            Enmp = fetmodel.Ep_nm_radial1d(p.ems, p.W1/2, int(n), int(m))
            gamma_nm = fetmodel.gamma_nm_NP(Enmp, p.alpha)
            alpha_nm = fetmodel.alpha_nm_NP(p.alpha,gamma_nm)
            Enm = fetmodel.E_nm_NP(p.alpha, gamma_nm)
            cur1 = func_FD0(EFs-Enm-e0, p.temp)
            cur2 = func_FD0(EFs-Enm-e0-Vds, p.temp)
            cur += cur1-cur2

    return (cur*2*const.elementary_charge/const.h*p.temp*const.Boltzmann)


####### functions
def func_FD0(ene, temp):
    """ Fermi-Dirac integral with order 0
    """
    ene0=ene*const.elementary_charge/(const.Boltzmann*temp)
    # print(ene0)
    if(ene0>100):
        return ene0
    else:
        return math.log(1+math.exp(ene0))

def determine_EFs_circ1d(p, Vds, Ids_cutoff):
    """
    function to find E_F at source in circular 1d NW
    """
    e0 = optimize.root_scalar(
        func_det_EFs_circ1d, args=(p, Vds, Ids_cutoff), x0=0.0, x1=0.15)
    if(e0.converged==True):
        return e0.root
    else:
        print("EFs convergence error !")
        print(e0)
        return 0


def func_det_EFs_circ1d(EFs, p, Vds, Ids_cutoff):
    """
    function to find E_F at source in circular 1d NW
    """
    ids0 = Ids_ballistic1d_circ1d(Vds, 0, p, EFs)/(math.pi*p.W1)
    return ids0-Ids_cutoff

def func_charge1d(Vgs, Vds,p):
    """
    1D charge density (assume E0 is obtained)
    """
    e0=get_E0_circ1d(p, Vgs, Vds)
    return(p.Ceff/const.elementary_charge*(e0+p.alpha_D*Vds+p.alpha_G*Vgs))


def determine_EFs_circ1dNP(p, Vds, Ids_cutoff):
    """
    function to find E_F at source in circular 1d NW
    """
    e0 = optimize.root_scalar(
        func_det_EFs_circ1dNP, args=(p, Vds, Ids_cutoff), x0=0.0, x1=0.15)
    if(e0.converged==True):
        return e0.root
    else:
        print("EFs convergence error !")
        print(e0)
        return 0


def func_det_EFs_circ1dNP(EFs, p, Vds, Ids_cutoff):
    """
    function to find E_F at source in circular 1d NW
    """
    ids0 = Ids_ballistic1d_circ1dNP(Vds, 0, p, EFs)/(math.pi*p.W1)
    return ids0-Ids_cutoff


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
    W1 = 30e-9 ## diameter
    W2 = W1
    alpha = fetmodel.alpha_NP(Eg, ems)
    Cox = fetmodel.Cox_radial(epsOX, tOX, W1/2)
    Cc = fetmodel.Cc_rect(epsS, W1, W1)
    # zI = radius-Delta,
    #   Delta = \int rho(r) r^2 dr / \int rho(r) r dr
    #    Cc = fetmodel.Cc_radial(epsS, zI, W1-zI)
    # See eq. in
    # IEEE TRANSACTIONS ON ELECTRON DEVICES, VOL. 55, NO. 1, JANUARY 2008 411
    # Modeling the Centroid and the Inversion Charge in
    # Cylindrical Surrounding Gate MOSFETs,
    # Including Quantum Effects
    # J. B. Roldán, Andrés Godoy, Francisco Gámiz, Senior Member, IEEE, and M. Balaguer
    nmax=5
    mmax=5
    print('Cox=', Cox,', Cc=', Cc)
    p=fetmodel.parameters_ballistic(EFermi=EFermi,
                                    alpha=alpha,
                                    Ceff=Cox*Cc/(Cox+Cc),
                                    ems=ems,
                                    W1=W1,
                                    W2=W2,
                                    nmax=nmax,
                                    mmax=mmax)
    p.output()

    Epnm=np.zeros(nmax*mmax)
    gamma_nm=np.zeros(nmax*mmax)
    alpha_nm=np.zeros(nmax*mmax)
    ems_nm=np.zeros(nmax*mmax)
    Enm=np.zeros(nmax*mmax)
    print('Energy Levels: n, m, parabollic, nonparabollic')
    i=0
    for n in np.arange(1, p.nmax+1):
        for m in np.arange(1,p.mmax+1):
            Epnm[i]=fetmodel.Ep_nm_radial1d(p.ems, p.W1/2, int(n), int(m))
            gamma_nm[i] = fetmodel.gamma_nm_NP(Epnm[i], alpha)
            alpha_nm[i] = fetmodel.alpha_nm_NP(alpha,gamma_nm[i])
            ems_nm[i] = fetmodel.ems_nm_NP(p.ems,gamma_nm[i])
            Enm[i] = fetmodel.E_nm_NP(alpha, gamma_nm[i])
            print(n,m,Epnm[i],Enm[i], gamma_nm[i],ems_nm[i],alpha_nm[i])
            i+=1

    plt.figure()
    left=min(Epnm)
    right=max(Epnm)
    plt.scatter(Epnm,Enm,label='energy')
    plt.legend(loc='best')
    
    left=min(Enm)
    right=max(Enm)
    plt.figure()
    plt.scatter(Enm,alpha_nm,label='alpha')
    plt.hlines([alpha], left,right, "blue", linestyles='dashed')
    plt.legend(loc='best')
    plt.figure()
    plt.scatter(Enm,ems_nm,label='ems')
    plt.hlines([ems], left, right, "blue", linestyles='dashed')
    plt.legend(loc='best')
    plt.show()

    EFs = 0.0
    Vgs = 0.1
    print(Vgs, Ids_ballistic1d_circ1d(0.5, Vgs, p, EFs)/(math.pi*p.W1), "uA/um")
    print(Vgs, Ids_ballistic1d_circ1dNP(0.5, Vgs, p, EFs)/(math.pi*p.W1), "uA/um")
    Vgs = 0
    print(Vgs, Ids_ballistic1d_circ1d(0.5, Vgs, p, EFs)/(math.pi*p.W1), "uA/um")
    print(Vgs, Ids_ballistic1d_circ1dNP(0.5, Vgs, p, EFs)/(math.pi*p.W1), "uA/um")
    Vgs = -0.1
    print(Vgs, Ids_ballistic1d_circ1d(0.5, Vgs, p, EFs)/(math.pi*p.W1), "uA/um")
    print(Vgs, Ids_ballistic1d_circ1dNP(0.5, Vgs, p, EFs)/(math.pi*p.W1), "uA/um")

    # # determine Vth/EFs so that the current at Vgs=0 is given as follows
    Vds = 0.5
    Ids_cutoff = 100e-9 * 1e6
    EFs = determine_EFs_circ1d(p, Vds, Ids_cutoff)
    EFs2 = determine_EFs_circ1dNP(p, Vds, Ids_cutoff)
    print(EFs, Ids_ballistic1d_circ1d(0.5, 0, p, EFs)/(math.pi*p.W1), "uA/um")
    print(EFs2, Ids_ballistic1d_circ1dNP(0.5, 0, p, EFs2)/(math.pi*p.W1), "uA/um")

    dVgs = 0.005
    Vgs = np.arange(0, 0.6, dVgs)
    Ids1 = np.empty_like(Vgs)
    Ids2 = np.empty_like(Vgs)
    for i, Vgs0 in enumerate(Vgs):
        Ids1[i] = Ids_ballistic1d_circ1d(Vds, Vgs0, p, EFs)/(math.pi*p.W1)*1e-3
        # Ids2[i] = Ids_ballistic1d_circ1dNP(Vds, Vgs0, p, EFs2)/(math.pi*p.W1)*1e-3

    gm1 = np.gradient(Ids1, dVgs)
    # gm2 = np.gradient(Ids2, dVgs)

    plt.plot(Vgs, Ids1, label='parabollic')
    # plt.plot(Vgs, Ids2, label='nonparabollic')
    plt.xlabel('Vgs - Vth (V)')
    plt.ylabel('Ids (mA/um)')
    plt.legend(loc='best')
    plt.show()

    plt.plot(Vgs, gm1, label='parabollic')
    # plt.plot(Vgs, gm2, label='nonparabollic')
    plt.xlabel('Vgs - Vth (V)')
    plt.ylabel('Transconductance (mS/um)')
    plt.legend(loc='best')
    plt.show()
