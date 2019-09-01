#!/usr/bin/env python
# -*- coding: utf-8 -*-

# test of density1d

import fetmodel
import numpy as np
import matplotlib.pyplot as plt

Eg = 0.36
ems = 0.0671
alpha = fetmodel.alpha_NP(Eg, ems)
temperature = 300

W1 = 10e-9
W2 = 8e-9

Enmp = fetmodel.Ep_nm_rect1d(ems, W1, W2, 1, 1)
gamma_nm = fetmodel.gamma_nm_NP(Enmp, alpha)
Enm = fetmodel.E_nm_NP(alpha, gamma_nm)
alpha_nm = fetmodel.alpha_nm_NP(alpha, gamma_nm)
ems_nm = fetmodel.ems_nm_NP(ems, gamma_nm)

print(alpha, gamma_nm, Enm, alpha_nm, ems_nm)

ene = np.arange(-0.1, 0.2, 0.005)
dens = np.empty_like(ene)
dens2 = np.empty_like(ene)
dens3 = np.empty_like(ene)

for i, e in enumerate(ene):
    dens[i] = fetmodel.density1d0(e, Enm, ems_nm, temperature)
    dens2[i] = fetmodel.density1d_NP0(
        e, Enm, alpha_nm, ems_nm, temperature)
    dens3[i] = fetmodel.density1d_rect1dNP_all0(e, alpha, ems, temperature,
                                             W1, W2, 2, 2)

plt.plot(ene, dens, ene, dens2, ene, dens3)
plt.show()
