#!/usr/bin/env python
# coding: utf-8

import fetmodel
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import scipy.constants as const


def density1d_circ1d_all0(EFermi, ems, temp, radius, nmax):
    """ electron denstiy in circular nanowire
    """
    n0 = 0
    nlist = np.arange(1, nmax+1, dtype=np.int64)
    for n in nlist:
        Enp = fetmodel.Ep_n_radial1d(ems, radius, int(n))
        nn = fetmodel.density1d0(EFermi, Enp, ems, temp)
        n0 += nn
    return(n0)


def func_e0_find_circ1d(E0, p, Vgs, Vds):
    """ function to find find energy of the top of the barrier in circular NW
    """
    n1d_S = density1d_circ1d_all0(p.EFermi - E0, p.ems, p.temp, p.W1, p.nmax)
    n1d_D = density1d_circ1d_all0(
        p.EFermi - E0 - Vds, p.ems, p.temp, p.W1, p.nmax)
    q0 = const.elementary_charge * (n1d_S + n1d_D) / (2 * p.Ceff)
    return E0 + (p.alpha_D * Vds + p.alpha_G * Vgs - q0)


def get_E0_circ1d(p, Vgs, Vds):
    """ get energy of the top of the barrier
    """
    e0 = optimize.root_scalar(
        func_e0_find_circ1d, args=(p, Vgs, Vds), x0=-0.1, x1=1)
    return e0.root


def func_FD0(ene, temp):
    """ Fermi-Dirac integral with order 0
    """
    return math.log(1+math.exp(ene*const.elementary_charge/(const.Boltzmann*temp)))


def func_current_circ1d(Vds, Vgs, p, EFs):
    """ current of ballistic FET with circular 1d NW
    """
    e0 = optimize.root_scalar(
        func_e0_find_circ1d, args=(p, Vgs, Vds), x0=-0.1, x1=1)
    # e0=get_E0_circ1d(p, Vgs, Vds)
    nlist = np.arange(1, p.nmax+1, dtype=np.int64)
    cur = 0
    for n in nlist:
        Enp = fetmodel.Ep_n_radial1d(p.ems, p.W1, int(n))
        cur1 = func_FD0(EFs-Enp-e0.root, p.temp)
        cur2 = func_FD0(EFs-Enp-e0.root-Vds, p.temp)
        cur += cur1-cur2

    return (cur*2*const.elementary_charge/6.63e-34*p.temp*const.Boltzmann)


def determine_EFs_circ1d(p, Vds, Ids):
    """ function to find E_F at source in circular 1d NW
    """
    e0 = optimize.root_scalar(
        func_det_EFs_circ1d, args=(p, Vds, Ids), x0=-0.1, x1=0.3)
    return e0.root


def func_det_EFs_circ1d(EFs, p, Vds, Ids):
    """ function to find E_F at source in circular 1d NW
    """
    ids0 = func_current_circ1d(Vds, 0, p, EFs)/(2*math.pi*p.W1)
    return Ids-ids0


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
    p.nmax = 4
    p.mmax = 2

    EFs = 0.095
    # Vgs = 0.1
    # print(Vgs, func_current_circ1d(0.5, Vgs, p, EFs)/(2*math.pi*p.W1), "uA/um")
    Vgs = 0
    # print(Vgs, func_current_circ1d(0.5, Vgs, p, EFs)/(2*math.pi*p.W1), "uA/um")
    # Vgs = -0.1
    # print(Vgs, func_current_circ1d(0.5, Vgs, p, EFs)/(2*math.pi*p.W1), "uA/um")

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
