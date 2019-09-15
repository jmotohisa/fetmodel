#!/usr/bin/env python
# coding: utf-8

import fetmodel
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import scipy.constants as const


def density1d_circ1d_all0(EFermi, ems, temp, radius, nmax, mmax):
    """ electron denstiy in circular nanowire
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


def func_e0_find_circ1d(E0, p, Vgs, Vds):
    """ function to find find energy of the top of the barrier in circular NW
    """
    n1d_S = density1d_circ1d_all0(p.EFermi - E0, p.ems, p.temp, p.W1, p.nmax,p.mmax)
    n1d_D = density1d_circ1d_all0(
        p.EFermi - E0 - Vds, p.ems, p.temp, p.W1, p.nmax,p.mmax)
    q0 = const.elementary_charge * (n1d_S + n1d_D) / (2 * p.Ceff)
    return E0 + (p.alpha_D * Vds + p.alpha_G * Vgs - q0)


def get_E0_circ1d(p, Vgs, Vds):
    """ get energy of the top of the barrier
    """
    e0 = optimize.root_scalar(
        func_e0_find_circ1d, args=(p, Vgs, Vds), x0=-0.0, x1=0.15)
    if(e0.converged==True):
        return e0.root
    else:
        print("E0 convergence error !")
        print(e0)
        return 0


def func_FD0(ene, temp):
    """ Fermi-Dirac integral with order 0
    """
    ene0=ene*const.elementary_charge/(const.Boltzmann*temp)
    # print(ene0)
    if(ene0>100):
        return ene0
    else:
        return math.log(1+math.exp(ene0))


def func_current_circ1d(Vds, Vgs, p, EFs):
    """ current of ballistic FET with circular 1d NW
    """
    # e00 = optimize.root_scalar(
    #     func_e0_find_circ1d, args=(p, Vgs, Vds), x0=-0.1, x1=1)
    # e0=e00.root
    e0=get_E0_circ1d(p, Vgs, Vds)
    nlist = np.arange(1, p.nmax+1, dtype=np.int64)
    mlist = np.arange(1, p.mmax+1, dtype=np.int64)
    cur = 0
    for m in mlist:
        for n in nlist:
            Enp = fetmodel.Ep_nm_radial1d(p.ems, p.W1, int(n), int(m))
            cur1 = func_FD0(EFs-Enp-e0, p.temp)
            cur2 = func_FD0(EFs-Enp-e0-Vds, p.temp)
            cur += cur1-cur2

    return (cur*2*const.elementary_charge/6.63e-34*p.temp*const.Boltzmann)


def determine_EFs_circ1d(p, Vds, Ids):
    """ function to find E_F at source in circular 1d NW
    """
    e0 = optimize.root_scalar(
        func_det_EFs_circ1d, args=(p, Vds, Ids), x0=0.0, x1=0.15)
    if(e0.converged==True):
        return e0.root
    else:
        print("EFs convergence error !")
        print(e0)
        return 0


def func_det_EFs_circ1d(EFs, p, Vds, Ids):
    """ function to find E_F at source in circular 1d NW
    """
    ids0 = func_current_circ1d(Vds, 0, p, EFs)/(2*math.pi*p.W1)
    return Ids-ids0

def func_charge1d(Vgs, Vds,p):
    """
    1D charge density (assume E0 is obtained)
    """
    e0=get_E0_circ1d(p, Vgs, Vds)
    return(p.Ceff/const.elementary_charge*(e0+p.alpha_D*Vds+p.alpha_G*Vgs))


if __name__ == '__main__':
    Eg = 0.36
    epsOX = 20
    epsS = 15.15
    tOX = 3e-9
    temperature = 300
    ems = 0.023
    W1 = 10e-9
    W2 = 8e-9
    alpha = fetmodel.alpha_NP(Eg, ems)
    Cox = fetmodel.Cox_radial(epsOX, tOX, W1)
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
    alpha_D = 0
    alpha_G = 1

    p = fetmodel.param_ballistic_new()
    p.ems = ems
    p.alpha = alpha
    p.W1 = W1
    p.W2 = W2

    p.EFermi = -0.1
    p.alpha_D = alpha_D
    p.alpha_G = alpha_G
    p.Ceff = Cox*Cc/(Cox+Cc)
    p.temp = temperature
    p.nmax = 2
    p.mmax = 4

    for n in np.arange(1, p.nmax+1):
        for m in np.arange(1,p.mmax+1):
            Epnm=fetmodel.Ep_nm_radial1d(p.ems, p.W1, int(n), int(m))
            print(n,m,Epnm)

    EFs = 0.095
    Vgs = 0.1
    print(Vgs, func_current_circ1d(0.5, Vgs, p, EFs)/(2*math.pi*p.W1), "uA/um")
    Vgs = 0
    print(Vgs, func_current_circ1d(0.5, Vgs, p, EFs)/(2*math.pi*p.W1), "uA/um")
    Vgs = -0.1
    print(Vgs, func_current_circ1d(0.5, Vgs, p, EFs)/(2*math.pi*p.W1), "uA/um")

    # # determine Vth/EFs so that the current at Vgs=0 is given as follows
    Vds = 0.5
    Ids_cutoff = 100e-9 * 1e6
    EFs = determine_EFs_circ1d(p, Vds, Ids_cutoff)
    print("EFs:", EFs)

    dVgs = 0.005
    Vgs = np.arange(0, 0.6, dVgs)
    Ids1 = np.empty_like(Vgs)
    for i, Vgs0 in enumerate(Vgs):
        Ids1[i] = func_current_circ1d(Vds, Vgs0, p, EFs)/(2*math.pi*p.W1)*1e-3

    gm1 = np.gradient(Ids1, dVgs)

    plt.plot(Vgs, Ids1, label='Ids1')
    plt.xlabel('Vgs - Vth (V)')
    plt.ylabel('Ids (mA/um)')
    plt.show()

    plt.plot(Vgs, gm1, label='gm')
    plt.xlabel('Vgs - Vth (V)')
    plt.ylabel('Transconductance (mS/um)')
    plt.show()
