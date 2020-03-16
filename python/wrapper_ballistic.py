#!/usr/bin/env python
# -*- coding: utf-8 -*-

import scipy.constants as const
import fetmodel
import numpy as np

def Ids_rect1dNP(Vds, Vgs, p, EFs):
    """
    p: class of parameters_ballistic
    (wrapper)
    """
    return(fetmodel.Ids_ballistic1d_rect1dNP0(Vds, Vgs, EFs, p.alpha_D, p.alpha_G, p.Ceff,
                                              p.alpha, p.ems, p.temp, p.W1, p.W2, int(p.nmax), int(p.mmax)))

def Ids_2d(Vds, Vgs, p, EFs):
    return fetmodel.Ids_ballistic2d0(Vds, Vgs, EFs,p.alpha_D, p.alpha_G,p.Ceff,p.ems, p.temp)
