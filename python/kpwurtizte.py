#!/usr/bin/env python
# -*- coding: utf-8 -*-

# unit: hbar=1 and m0=1

import math
import numpy as np
import matplotlib.pyplot as plt
# import ballistic
from scipy import optimize


def func_P1_2(mz, Eg, Delta1, Delta2, Delta3):
    return (1/mz-1)*((Eg+Delta1+Delta2)*(Eg+2*Delta2)-2*Delta3**2)/(Eg+2*Delta2)/2


def func_P2_2(mt, Eg, Delta1, Delta2, Delta3):
    return (1/mt-1)*(Eg*((Eg+Delta1+Delta2)*(Eg+2*Delta2)-2*Delta3**2) / ((Eg+Delta1+Delta2)*(Eg+Delta2)-Delta3**2))/2


def Eck(kx, ky, kz, mz, mt):
    return ((kx**2+ky**2)/(2 * mt) + kz**2/(2 * mz))


def eigenWZ(kx, ky, kz, P12, P22, Eg, Delta1, Delta2, Delta3):
    c0 = kz**2*P12*Delta1**2 - kz**2*P12 * Delta2**2 - (kx**2 + ky**2)*P22*Delta3**2 + 2*Eg*Delta1*Delta3**2 + \
        2*Delta1**2*Delta3**2 + 2*Eg*Delta2*Delta3**2 + 4 * \
        Delta1*Delta2*Delta3**2 + 2*Delta2**2*Delta3**2
    c1 = - 2*kz**2*P12*Delta1 + (-kx**2 - ky**2)*P22*Delta1 + Eg*Delta1**2 + Delta1**3 + Delta1**2*Delta2 - \
        Eg*Delta2**2 - Delta1*Delta2**2 - Delta2**3 - 2*Eg * \
        Delta3**2 - 4*Delta1*Delta3**2 - 4*Delta2*Delta3**2
    c2 = (kz**2*P12 + (kx**2 + ky**2)*P22 - 2*Eg*Delta1 - 3 *
          Delta1**2 - 2*Delta1*Delta2 + Delta2**2 + 2*Delta3**2)
    c3 = (Eg + 3*Delta1 + Delta2)
    c4 = -1
    return np.roots([c4, c3, c2, c1, c0])


def eigenZB(Eg, deltaSO, kx, ky, kz, P02):
    deltaSO3 = deltaSO/3
    c0 = deltaSO3**2*(-P02*k**2 + 2*Eg*deltaSO3)
    c1 = -(deltaSO3**2*(3*Eg + 2*deltaSO3))
    c2 = k**2*P02 + 3*deltaSO3**2
    c3 = Eg
    c4 = -1
    return np.roots([c4, c3, c2, c1, c0])


# GaN
Eg = 3.4
Delta_cr = 0.017
Delta_so = 0.01
mt = 0.2
mz = 0.2

Delta1 = Delta_cr
Delta2 = Delta_so/3
Delta3 = Delta_so/3
P12 = func_P1_2(mz, Eg, Delta1, Delta2, Delta3)
P22 = func_P2_2(mt, Eg, Delta1, Delta2, Delta3)
print(math.sqrt(P12), math.sqrt(P22))

print(eigenWZ(0, 0, 0, P12, P22, Eg, Delta1, Delta2, Delta3))

kz = np.arange(0, 0.2, 0.001)
e1 = np.empty_like(kz)
e2 = np.empty_like(kz)
for i, kz0 in enumerate(kz):
    # eig = [0, 0, 0, 0]
    eig = eigenWZ(0, 0, kz0, P12, P22, Eg, Delta1, Delta2, Delta3)
    e1[i] = eig[0]
    e2[i] = Eck(0, 0, kz0, mz, mt) + Eg + Delta1 + Delta2

plt.plot(kz, e1, label='nonparab')
plt.plot(kz, e2, label='parabollic')
plt.legend(loc='best')
plt.show()
