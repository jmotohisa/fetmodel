#!/usr/bin/env python
# -*- coding: utf-8 -*-

# test of ballisti 2d FET

import fetmodel
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import ballistic2d

Eg = 0.36
epsOX = 8.5
epsS = 8.9
tOX = 20e-9
temperature = 300
ems = 0.067
W1 = 10e-9
W2 = 8e-9
alpha = fetmodel.alpha_NP(Eg, ems)
Cox = epsOX*8.85e-12/tOX
Cc = math.sqrt(1.6e-19*1e21/(2*epsS*8.86e-12*1))
alpha_D = 0
alpha_G = 1

p = fetmodel.param_ballistic_new()
p.ems = ems
p.alpha = alpha
p.W1 = W1
p.W2 = W2

p.EFermi = -0.1
p.VDS = 0
p.VGS = 0
p.alpha_D = alpha_D
p.alpha_G = alpha_G
p.Ceff = Cox*Cc/(Cox+Cc)
p.temp = temperature
p.nmax = 2
p.mmax = 2

# Vgs = 0
# Vds = 0
# e = np.arange(-0.1, 0.1, 0.005)
# e0 = np.empty_like(e)
# for i, e00 in enumerate(e):
#     e0[i] = ballistic2d.func_e0_find(e00, p, Vgs, Vds)

# plt.plot(e, e0)
# plt.xlabel('Energy (eV)')
# plt.ylabel('func_E0_find')
# plt.show()

# Vds = 0
# Vgs = np.arange(0, 1, 0.01)
# e0 = np.empty_like(Vgs)

# for i, Vgs0 in enumerate(Vgs):
#     e0[i] = ballistic2d.get_E0(p, Vgs0, Vds)

# plt.plot(Vgs, e0)
# plt.xlabel('Gate Voltage (V)')
# plt.ylabel('Top of the barrier height (eV)')
# plt.show()

# Vds = np.arange(0, 1, 0.01)
# Ids1 = np.empty_like(Vds)
# Ids2 = np.empty_like(Vds)

# for i, Vds0 in enumerate(Vds):
#     Ids1[i] = ballistic2d.Ids_ballistic2d(Vds0, -0.1, p, 0)
#     Ids2[i] = ballistic2d.Ids_ballistic2d(Vds0, 0, p, 0)

# plt.plot(Vds, Ids1, label='Ids1')
# plt.plot(Vds, Ids2, label='Ids2')
# plt.xlabel('Drain Voltage (V)')
# plt.ylabel('Drain Current (mA/um)')
# plt.legend(loc='best')
# plt.show()

Vgs = np.arange(-0.1, 1, 0.01)
Ids1 = np.empty_like(Vgs)
Ids2 = np.empty_like(Vgs)
Vds = 0.5

for i, Vgs0 in enumerate(Vgs):
    Ids1[i] = ballistic2d.Ids_ballistic2d(Vds, Vgs0, p, 0)*1e-3
    Ids2[i] = ballistic2d.Ids_ballistic2d0(Vds, Vgs0, p, 0)*1e-3

fig, ax = plt.subplots()
ax.plot(Vgs, Ids1, label='Ids1')
ax.plot(Vgs, Ids2, label='Ids2')
ax.set_yscale("log")
plt.xlabel('Gate Voltage (V)')
plt.ylabel('Drain Current (mA/um)')
plt.legend(loc='best')
plt.show()

# for i, ids0 in enumerate(Ids1):
#     print(ids0/Ids2[i])
