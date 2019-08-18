#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import numpy as np
import matplotlib.pyplot as plt
import ballistic
from scipy import optimize

Eg = 0.4
ems = 0.067
alpha = ballistic.alphaNP00(Eg, ems)
temperature = 300
epsOX = 8.5
epsS = 8.9
W1 = 10e-9
W2 = 8e-9
tOX = 20e-9
Cox = ballistic.Cox_rect(epsOX, tOX, W1, W2)
Cc = ballistic.Cc_rect(epsS, W1, W2)
p = ballistic.param_ballistic_new()

p.EFermi = 0
p.VDS = 0
p.VGS = 0
p.alpha_D = 0
p.alpha_G = 1
p.Ceff = Cox*Cc/(Cox+Cc)
p.alpha = alpha
p.ems = ems
p.temp = temperature
p.W1 = W1
p.W2 = W2
p.nmax = 2
p.mmax = 2

alpha_D = 0
alpha_G = 1

# print(ballistic.find_E0(p))


# ballistic.density1d_all00(e, alpha, ems, temperature,
#                           W1, W2, 2, 2)


def func_e0_find(E0, EFermi, Vds, Vgs, alpha_D, alpha_G, Ceff):
    n1d_S = ballistic.density1d_all00(
        E0-EFermi, alpha, ems, temperature, W1, W2, 2, 2)
    n1d_D = ballistic.density1d_all00(
        E0-EFermi-Vds, alpha, ems, temperature, W1, W2, 2, 2)
    q0 = 1.6e-19*(n1d_S+n1d_D)/(2*Ceff)
    return E0-(alpha_D*Vds+alpha_G*Vgs-q0)


EFermi = 0
Vds = 0
Vgs = 0


e = np.arange(-0.1, 0.1, 0.005)
e0 = np.empty_like(e)
for i, e00 in enumerate(e):
    e0[i] = func_e0_find(e00, EFermi, Vds, Vgs, alpha_D, alpha_G, p.Ceff)

plt.plot(e, e0)
plt.show()

Vgs = np.arange(0, 1, 0.01)
e0 = np.empty_like(Vgs)

for i, Vgs0 in enumerate(Vgs):
    e00 = optimize.root_scalar(func_e0_find, args=(
        p.EFermi, Vds, Vgs0, alpha_D, alpha_G, p.Ceff), x0=-0.1, x1=1)
    e0[i] = e00.root

plt.plot(Vgs, e0)
plt.show()


def func_FD0(ene, temp):
    return math.log(1+math.exp(ene*1.6e-19/(1.38e-23*temp)))


def func_current1D(Vgs, Vds, p):
    e0 = optimize.root_scalar(func_e0_find, args=(
        EFermi, Vds, Vgs, p.alpha_D, p.alpha_G, p.Ceff), x0=-0.1, x1=1)
    nlist = np.arange(1, p.nmax+1, dtype=np.int64)
    mlist = np.arange(1, p.mmax+1, dtype=np.int64)
    cur = 0
    for n in nlist:
        for m in mlist:
            gamma_nm = ballistic.gamma_nm00(
                alpha, ems, p.W1, p.W2, int(n), int(m))
            Enm = ballistic.E_nm0(alpha, gamma_nm)
            cur1 = func_FD0(EFermi-Enm-e0.root, temperature)
            cur2 = func_FD0(EFermi-Enm-e0.root-Vds, temperature)
            cur += cur1-cur2

    return (cur*2*1.6e-19/6.63e-34)


Vds = np.arange(0, 1, 0.01)
Ids1 = np.empty_like(Vds)
Ids2 = np.empty_like(Vds)

for i, Vds0 in enumerate(Vds):
    Ids1[i] = func_current1D(-0.1, Vds0, p)
    Ids2[i] = func_current1D(0, Vds0, p)

plt.plot(Vds, Ids1, Vds, Ids2)
plt.show()
