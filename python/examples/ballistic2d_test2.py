#!/usr/bin/env python
# -*- coding: utf-8 -*-

# test of ballisti 2d FET

import pyfet
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

Eg = 0.36
epsOX = 8.5
epsS = 8.9
tOX = 20e-9
temperature = 300
ems = 0.067
W1 = 10e-9
W2 = 8e-9
alpha = pyfet.alpha_NP(Eg, ems)
Cox = epsOX*8.85e-12/tOX
Cc = math.sqrt(epsS)
alpha_D = 0
alpha_G = 1

p = pyfet.param_ballistic_new()
p.ems = ems
# p.alpha = alpha
# p.W1 = W1
# p.W2 = W2

p.EFermi = -0.1
p.alpha_D = alpha_D
p.alpha_G = alpha_G
p.Ceff = Cox*Cc/(Cox+Cc)
p.temp = temperature
# p.nmax = 2
# p.mmax = 2


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


Vds = np.arange(0, 1, 0.01)
Ids1 = np.empty_like(Vds)
Ids2 = np.empty_like(Vds)

for i, Vds0 in enumerate(Vds):
    Ids1[i] = pyfet.Ids_ballistic2d(Vds0, -0.1, p, 0)
    Ids2[i] = pyfet.Ids_ballistic2d(Vds0, 0, p, 0)

plt.plot(Vds, Ids1, label='Ids1')
plt.plot(Vds, Ids2, label='Ids2')
plt.xlabel('Drain Voltage (V)')
plt.ylabel('Drain Current (A/m)')
plt.legend(loc='best')
plt.show()
