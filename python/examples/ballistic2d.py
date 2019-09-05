#!/usr/bin/env python
# -*- coding: utf-8 -*-

import fetmodel
import math
import fdint
import scipy.constants as const
from scipy import optimize


def func_e0_find(E0, p, Vgs, Vds):
    n2d_S = fetmodel.density2d0(p.EFermi - E0, 0, p.ems, p.temp)
    n2d_D = fetmodel.density2d0(p.EFermi - E0 - Vds, 0, p.ems, p.temp)
    q0 = const.elementary_charge * (n2d_S + n2d_D) / (2 * p.Ceff)
    return E0 + (p.alpha_D * Vds + p.alpha_G * Vgs - q0)


def get_E0(p, Vgs, Vds):
    e0 = optimize.root_scalar(func_e0_find, args=(p, Vgs, Vds), x0=-1, x1=1)
    return e0.root


def Ids_ballistic2d(Vds, Vgs, p, EFs):
    ee0 = get_E0(p, Vgs, Vds)
    return(fetmodel.Ids_ballistic2d0_E0(ee0, Vds, Vgs, EFs,
                                        p.EFermi, p.alpha_D, p.alpha_G, p.Ceff, p.ems, p.temp))


def Ids_ballistic2d0(Vds, Vgs, p, EFs):
    ee0 = get_E0(p, Vgs, Vds)
    n2d = p.ems*const.electron_mass / \
        (math.pi*const.hbar**2)*p.temp*const.Boltzmann
    beta = const.elementary_charge/(p.temp*const.Boltzmann)
    ns = (fdint.fdk(0, beta*(EFs-ee0))+fdint.fdk(0, beta*(EFs-ee0-Vds)))/2*n2d
    v0 = math.sqrt(2*p.temp*const.Boltzmann /
                   (math.pi*p.ems*const.electron_mass))
    v1 = fdint.fdk(0.5, beta*(EFs-ee0))*(2/math.sqrt(math.pi))
    v2 = fdint.fdk(0, beta*(EFs-ee0))
    vinj = v0*v1/v2
    f1 = 1 - fdint.fdk(0.5, beta*(EFs-ee0-Vds))/fdint.fdk(0.5, beta*(EFs-ee0))
    f2 = 1 + fdint.fdk(0, beta*(EFs-ee0-Vds))/fdint.fdk(0, beta*(EFs-ee0))
    return const.elementary_charge*ns*vinj*f1/f2
