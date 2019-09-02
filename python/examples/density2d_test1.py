#!/usr/bin/env python
# -*- coding: utf-8 -*-

# test of density2d

import fetmodel
import numpy as np
import math
import matplotlib.pyplot as plt

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

p.EFermi = 0
p.alpha_D = alpha_D
p.alpha_G = alpha_G
p.Ceff = Cox*Cc/(Cox+Cc)
p.temp = temperature

EFermi_list = np.arange(-0.15, 0.1, 0.005)
ns = np.empty_like(EFermi_list)
for i, EF0 in enumerate(EFermi_list):
    ns[i] = fetmodel.density2d0(EF0, 0, p.ems, p.temp)

plt.plot(EFermi_list, ns)
plt.xlabel('Fermi Energy (eV)')
plt.ylabel('2DEG density (m^-2)')
plt.show()
