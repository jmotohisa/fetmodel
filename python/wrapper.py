#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Wrapper functions for fetmodels
"""

import scipy.constants as const
import fetmodel
import numpy as np
import matplotlib.pyplot as plt

## from denstity1d.c
def density1d_rect1d_all(EFermi,p):
    """
    p: class of parameters_ballistic
    Wrapper function for density1d_rect1d_all0
    """
    return fetmodel.density1d_rect1d_all0(Efermi, p.ems, p.temp, p.W1, p.W2, int(p.nmax), int(p.mmax))

def density1d_rect1dNP_all(EFermi, p):
    """
    p: class of parameters_ballistic
    Wrapper function for density1d_rect1dNP_all0
    """
    return fetmodel.density1d_rect1dNP_all0(EFermi, p.alpha, p.ems, p.temp, p.W1, p.W2, int(p.nmax), int(p.mmax))

## from ballistic.c
## 1D

def func_for_findroot_E0_rect1dNP(ene0, Vds, Vgs, EFs, p):
    """
    C-implementation of the function to find top of the barrier (wrapper)
    p: class of parameters_ballistic
    Nonparabolic band
    """
    return fetmodel.func_for_findroot_E0_rect1dNP0(ene0, EFs, Vds, Vgs, p.alpha_D, p.alpha_G,
                                                   p.Ceff, p.alpha, p.ems,
                                                   p.temp, p.W1, p.W2, int(p.nmax), int(p.mmax))


def E0_rect1dNP_root(Vds,Vgs,EFermi,p):
    """
    Get top of the barrier (C-implementation, wrapper)
    p: class of parameters_ballistic
    Nonparabolic band
    """
    return fetmodel.E0_rect1dNP_root0(EFermi, Vds, Vgs, p.alpha_D, p.alpha_G, p.Ceff, p.alpha,
                                      p.ems, p.temp,p.W1,p.W2,int(p.nmax),int(p.mmax))


def Ids_ballistic1d_rect1dNP(Vds, Vgs, p, EFs):
    """
    p: class of parameters_ballistic
    (wrapper)
    """
    return(fetmodel.Ids_ballistic1d_rect1dNP0(Vds, Vgs, EFs, p.alpha_D, p.alpha_G, p.Ceff,
                                              p.alpha, p.ems, p.temp, p.W1, p.W2, int(p.nmax), int(p.mmax)))
    
## 2D

def func_for_findroot_E0_2d(ene0, Vds, Vgs, EFs, p):
    """
    C-implementation of the function to find top of the barrier in 2D ballistic FET (wrapper) 
    p: class of parameters_ballistic
    Single paraolic band
    """
    return fetmodel.func_for_findroot_E0_2d0(ene0, EFs, Vds, Vgs,
                                             p.alpha_D, p.alpha_G, p.Ceff,  p.ems, p.temp)

def E0_2d_root(Vds, Vgs, EFermi, p):
    """
    Get top of the barrier in 2D ballistic FET (C-implementation, wrapper)
    p: class of parameters_ballistic
    Single paraolic band
    """
    return fetmodel.E0_2d_root0(EFermi,Vds, Vgs, p.alpha_D, p.alpha_G, p.Ceff,p.ems, p.temp)
                
def Ids_ballistic2d_E0(E0,Vds, Vgs, p, EFs):
    return fetmodel.Ids_ballistic2d0_E0(E0, Vds, Vgs, EFs,p.alpha_D, p.alpha_G,Ceff,ems, p.temp)

def Ids_ballistic2d(Vds, Vgs, p, EFs):
    return fetmodel.Ids_ballistic2d0(Vds, Vgs, EFs,p.alpha_D, p.alpha_G,p.Ceff,p.ems, p.temp)


## from mos1d
def func_Qcharge1dMOS0(psis, ni, eps_semi, Na, temp):
    pp0 = Na
    np0 = ni**2/Na
    return fetmodel.func_QsMOS1D(psis,np0,pp0,eps_semi,temp)

## from cMOSFET
def Ids_cMOSFET(Vds, Vgs, p, ps=fetmodel.param_solver()):
    """
    p: class of parameters_cMOSFET
    """
    return(fetmodel.func_Ids_cMOSFET(Vds,Vgs,p,ps))


def Ids2_cMOSFET(Vds, Vgs, p):
    """
    p: class of parameters_cMOSFET
    """
    return(fetmodel.func_Ids2_cMOSFET(Vds,Vgs,p))
