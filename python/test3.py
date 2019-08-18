#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import numpy as np
import matplotlib.pyplot as plt
import ballistic
from scipy import optimize

Eg = 0.036
epsOX = 8.5
epsS = 8.9
tOX = 20e-9

p = ballistic.param_ballistic_new()
p.ems = 0.067
p.alpha = ballistic.alphaNP00(Eg, p.ems)
p.W1 = 10e-9
p.W2 = 8e-9

Cox = ballistic.Cox_rect(epsOX, tOX, p.W1, p.W2)
Cc = ballistic.Cc_rect(epsS, p.W1, p.W2)

p.EFermi = -0.1
# p.VDS = 0
# p.VGS = 0
p.alpha_D = 0
p.alpha_G = 1
p.Ceff = Cox*Cc/(Cox+Cc)
p.temp = 300
p.nmax = 2
p.mmax = 2


def func_e0_find(E0, p, Vgs, Vds):
    n1d_S = ballistic.density1d_all00(
        p.EFermi - E0, p.alpha, p.ems, p.temp, p.W1, p.W2, p.nmax, p.mmax)
    n1d_D = ballistic.density1d_all00(
        p.EFermi - E0 - Vds, p.alpha, p.ems, p.temp, p.W1, p.W2, p.nmax, p.mmax)
    q0 = 1.6e-19 * (n1d_S + n1d_D) / (2 * p.Ceff)
    return E0 + (p.alpha_D * Vds + p.alpha_G * Vgs - q0)


def get_E0(p, Vgs, Vds):
    e0 = optimize.root_scalar(func_e0_find, args=(p, Vgs, Vds), x0=-0.1, x1=1)
    return e0.root


# Vgs = 0
# Vds = 1
# e = np.arange(-1, 1, 0.005)
# e0 = np.empty_like(e)
# for i, e00 in enumerate(e):
#     e0[i] = func_e0_find(e00, p, Vgs, Vds)

# plt.plot(e, e0)
# plt.show()

# Vds = 0
# Vgs = np.arange(-0.1, 1, 0.01)
# e0 = np.empty_like(Vgs)

# for i, Vgs0 in enumerate(Vgs):
#     e0[i] = get_E0(p, Vgs0, Vds)

# plt.plot(Vgs, e0)
# plt.show()


def func_FD0(ene, temp):
    return math.log(1+math.exp(ene*1.6e-19/(1.38e-23*temp)))


def func_current1D(Vgs, Vds, p, EFs):
    e0 = get_E0(p, Vgs, Vds)
    nlist = np.arange(1, p.nmax+1, dtype=np.int64)
    mlist = np.arange(1, p.mmax+1, dtype=np.int64)
    cur = 0
    for n in nlist:
        for m in mlist:
            gamma_nm = ballistic.gamma_nm00(
                p.alpha, p.ems, p.W1, p.W2, int(n), int(m))
            Enm = ballistic.E_nm0(p.alpha, gamma_nm)
            cur1 = func_FD0(EFs-Enm-e0, p.temp)
            cur2 = func_FD0(EFs-Enm-e0-Vds, p.temp)
            cur += cur1-cur2

    return (cur*2*1.6e-19/6.63e-34*p.temp*1.38e-23)


Vds = np.arange(0, 1, 0.01)
Ids1 = np.empty_like(Vds)
Ids2 = np.empty_like(Vds)

for i, Vds0 in enumerate(Vds):
    Ids1[i] = func_current1D(-0.1, Vds0, p, 0)/(2*(p.W1+p.W2))
    Ids2[i] = func_current1D(0, Vds0, p, 0)/(2*(p.W1+p.W2))

plt.plot(Vds, Ids1, label='Ids1')
plt.plot(Vds, Ids2, label='Ids2')
plt.xlabel('Drain Voltage (V)')
plt.ylabel('Drain Current (A/m)')
plt.legend(loc='best')
plt.show()
