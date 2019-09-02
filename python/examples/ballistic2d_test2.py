#!/usr/bin/env python
# -*- coding: utf-8 -*-

# test of ballisti 2d FET

import fetmodel
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
p.alpha_D = alpha_D
p.alpha_G = alpha_G
p.Ceff = Cox*Cc/(Cox+Cc)
p.temp = temperature
p.nmax = 2
p.mmax = 2


Vds = np.arange(0, 1, 0.01)
Ids1 = np.empty_like(Vds)
Ids2 = np.empty_like(Vds)

for i, Vds0 in enumerate(Vds):
    Ids1[i] = fetmodel.Ids_ballistic2d0(
        Vds0, -0.1, 0, p.EFermi, p.alpha_D, p.alpha_G, p.Ceff, p.ems, p.temp)*1e-3
    Ids2[i] = fetmodel.Ids_ballistic2d0(
        Vds0, 0, 0, p.EFermi, p.alpha_D, p.alpha_G, p.Ceff, p.ems, p.temp)*1e-3

plt.plot(Vds, Ids1, label='Ids1')
plt.plot(Vds, Ids2, label='Ids2')
plt.xlabel('Drain Voltage (V)')
plt.ylabel('Drain Current (mA/um)')
plt.legend(loc='best')
plt.show()
