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
def density1d_rect1d_all(p):
    """
    p: class of parameters_ballistic
    Wrapper function for density1d_rect1d_all0
    """
    return fetmodel.density1d_rect1d_all0(p.EFermi, p.ems, p.temp, p.W1, p.W2, int(p.nmax), int(p.mmax))

def density1d_rect1dNP_all(p):
    """
    p: class of parameters_ballistic
    Wrapper function for density1d_rect1dNP_all0
    """
    return fetmodel.density1d_rect1dNP_all0(p.EFermi, p.alpha, p.ems, p.temp, p.W1, p.W2, int(p.nmax), int(p.mmax))

## from ballistic.c
## 1D

def func_for_findroot_E0_rect1dNP(ene0, Vds, Vgs, p):
    """
    C-implementation of the function to find top of the barrier (wrapper)
    p: class of parameters_ballistic
    Nonparabolic band
    """
    return fetmodel.func_for_findroot_E0_rect1dNP0(ene0, p.EFermi, Vds, Vgs, p.alpha_D, p.alpha_G,
                                                   p.Ceff, p.alpha, p.ems,
                                                   p.temp, p.W1, p.W2, int(p.nmax), int(p.mmax))


def E0_rect1dNP_root(Vds,Vgs,p):
    """
    Get top of the barrier (C-implementation, wrapper)
    p: class of parameters_ballistic
    Nonparabolic band
    """
    return fetmodel.E0_rect1dNP_root0(p.EFermi, Vds, Vgs, p.alpha_D, p.alpha_G, p.Ceff, p.alpha,
                                      p.ems, p.temp,p.W1,p.W2,int(p.nmax),int(p.mmax))


def Ids_ballistic1d_rect1dNP(Vds, Vgs, p, EFs):
    """
    p: class of parameters_ballistic
    (wrapper)
    """
    return(fetmodel.Ids_ballistic1d_rect1dNP0(Vds, Vgs, EFs, p.EFermi, p.alpha_D, p.alpha_G, p.Ceff,
                                              p.alpha, p.ems, p.temp, p.W1, p.W2, int(p.nmax), int(p.mmax)))
    

