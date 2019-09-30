#!/usr/bin/env python
# -*- coding: utf-8 -*-

import fetmodel
import numpy as np
import matplotlib.pyplot as plt
import math
import fdint
import scipy.constants as const
from scipy import optimize


def func_for_findroot_E0_2d(ene0, Vds,Vgs, EFs, p):
    """
    Python implementation of the function to find top of the barrier in 2D ballistic FET
    Parabolic band
    """
    n2d_S = fetmodel.density2d0(EFs - ene0, 0, p.ems, p.temp)
    n2d_D = fetmodel.density2d0(EFs - ene0 - Vds, 0, p.ems, p.temp)
    q0 = const.elementary_charge * (n2d_S + n2d_D) / (2 * p.Ceff)
    return ene0 + (p.alpha_D * Vds + p.alpha_G * Vgs - q0)


def E0_2d_root(Vds, Vgs, EFs, p, left=-0.5, right=0.5):
    """
    Get top of the barrier in 2D ballitic FET (Python implementation)
    Parabolic band
    """
    e0 = optimize.root_scalar(func_for_findroot_E0_2d, args=(Vds,Vgs, EFs, p), x0=left, x1=right)
    if(e0.converged==True):
        return e0.root
    else:
        print("EFs convergence error !")
        print(e0)
        return 0

# def get_E0_approx(p,Vgs,Vds):
#     c=const.elementary_charge**2*p.ems*const.electron_mass/(p.Ceff*math.pi*const.hbar**2)
#     return -(p.alpha_G*Vgs+p.alpha_D*Vds)/(1+c) + c/(1+c)*(p.EFermi-Vds/2)

def get_E0_approx(Vds, Vgs, EFs, P):
    n2dos = p.ems*const.electron_mass / \
            (math.pi*const.hbar**2)*p.temp*const.Boltzmann
    e0=-p.alpha_G*Vgs-p.alpha_D*Vds-const.elementary_charge*n2dos/p.Ceff
    return e0
    
    
def Ids_ballistic2d(Vds, Vgs, p, EFs):
    ene0 = E0_2d_root(Vds, Vgs, EFs, p)
    return(fetmodel.Ids_ballistic2d0_E0(ene0, Vds, Vgs, EFs,
                                        p.alpha_D, p.alpha_G, p.Ceff, p.ems, p.temp))


def Ids_ballistic2d_0(Vds, Vgs, p, EFs):
    ene0 = E0_2d_root(Vds, Vgs, EFs, p)
    n2d = p.ems*const.electron_mass / \
        (math.pi*const.hbar**2)*p.temp*const.Boltzmann
    beta = const.elementary_charge/(p.temp*const.Boltzmann)
    ns = (fdint.fdk(0, beta*(EFs-ene0))+fdint.fdk(0, beta*(EFs-ene0-Vds)))/2*n2d
    v0 = math.sqrt(2*p.temp*const.Boltzmann /
                   (math.pi*p.ems*const.electron_mass))
    v1 = fdint.fdk(0.5, beta*(EFs-ene0))*(2/math.sqrt(math.pi))
    v2 = fdint.fdk(0, beta*(EFs-ene0))
    vinj = v0*v1/v2
    f1 = 1 - fdint.fdk(0.5, beta*(EFs-ene0-Vds))/fdint.fdk(0.5, beta*(EFs-ene0))
    f2 = 1 + fdint.fdk(0, beta*(EFs-ene0-Vds))/fdint.fdk(0, beta*(EFs-ene0))
    return const.elementary_charge*ns*vinj*f1/f2

def determine_EFs_ballistic2d(p, Vds, Ids_cutoff):
    """ function to find E_F at source in Ballistic 2D FET
    """
    e0 = optimize.root_scalar(
        func_det_EFs_ballistic2d, args=(p, Vds, Ids_cutoff), x0=0.0, x1=0.15)
    if(e0.converged==True):
        return e0.root
    else:
        print("EFs convergence error !")
        print(e0)
        return 0

def func_det_EFs_ballistic2d(EFs, p,Vds,Ids_cutoff):
    ids0=fetmodel.Ids_ballistic2d(Vds, 0, p, EFs)
    return Ids_cutoff-ids0


def n2D(Vds, Vgs, p, EFs):
    ene0 = E0_2d_root(Vds, Vgs, EFs, p)
    n2dos = p.ems*const.electron_mass / \
            (math.pi*const.hbar**2)*p.temp*const.Boltzmann
    beta = const.elementary_charge/(p.temp*const.Boltzmann)
    ns = (fdint.fdk(0, beta*(EFs-ene0))+fdint.fdk(0, beta*(EFs-ene0-Vds)))/2*n2dos
    return ns

def n2D_approx(Vds,Vgs,p,EFs):
    ene0 = E0_2d_root(Vds, Vgs, EFs, p)
    n2dos = p.ems*const.electron_mass / \
            (math.pi*const.hbar**2)
    ns = n2dos*const.elementary_charge*(EFs-ene0-Vds/2)
    if(ns>0):
        return ns
    else:
        return 0


if __name__ == '__main__':
    Eg = 0.36
    epsOX = 8.5
    epsS = 8.9
    tOX = 20e-9
    temperature = 300
    ems = 0.067
    # W1 = 10e-9
    # W2 = 8e-9
    # alpha = fetmodel.alpha_NP(Eg, ems)
    Cox = epsOX*8.85e-12/tOX
    Cc = math.sqrt(1.6e-19*1e21/(2*epsS*8.86e-12*1))
    # Ceff=Cox*Cc/(Cox+Cc);
    Ceff=Cox
    # alpha_D = 0
    # alpha_G = 1
    print(Cox,Cc)
    p=fetmodel.parameters_ballistic(Ceff=Ceff,
                                    ems=ems,
                                    )
    p.output()
    print("Test of ballistic2d: parabolic band")

    Vgs = np.arange(-0.1, 1, 0.01)
    Ids1 = np.empty_like(Vgs)
    Ids2 = np.empty_like(Vgs)
    Vds = 0.5

    for i, Vgs0 in enumerate(Vgs):
        Ids1[i] = fetmodel.Ids_ballistic2d(Vds, Vgs0, p, 0)*1e-3
        Ids2[i] = Ids_ballistic2d(Vds, Vgs0, p, 0)*1e-3

    fig, ax = plt.subplots()
    ax.plot(Vgs, Ids1, label='Ids1')
    ax.plot(Vgs, Ids2, label='Ids2')
    ax.set_yscale("log")
    plt.xlabel('Gate Voltage (V)')
    plt.ylabel('Drain Current (mA/um)')
    plt.legend(loc='best')
    plt.show()

    Ids_cutoff = 100e-9*1e6
    Vgs0=0
    print("Determine EFs using fetmodel")
    EFs = determine_EFs_ballistic2d(p, Vds, Ids_cutoff)
    print("Source Fermi Energy for appropriate Vth (eV):", EFs)
    Ids0 = fetmodel.Ids_ballistic2d(Vds, Vgs0, p, EFs)
    print("Actual Ids at Vgs=0 (nA/um): ", Ids0*1e3)
    print("")

    Vds = np.arange(0, 1, 0.01)
    Ids1 = np.empty_like(Vds)
    Ids2 = np.empty_like(Vds)
    
    for i, Vds0 in enumerate(Vds):
        Ids1[i] = fetmodel.Ids_ballistic2d(Vds0, 0.5, p, EFs)*1e-3
        Ids2[i] = fetmodel.Ids_ballistic2d(Vds0, 0, p, EFs)*1e-3

    plt.plot(Vds, Ids1, label='Ids1')
    plt.plot(Vds, Ids2, label='Ids2')
    plt.xlabel('Drain Voltage (V)')
    plt.ylabel('Drain Current (mA/um)')
    plt.legend(loc='best')
    plt.show()
