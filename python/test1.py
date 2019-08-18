#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import ballistic

Eg = 0.4
ems = 0.067
alpha = ballistic.alphaNP00(Eg, ems)
temperature = 300

W1 = 10e-9
W2 = 8e-9

gamma_nm = ballistic.gamma_nm00(alpha, ems, W1, W2, 1, 1)
Enm = ballistic.E_nm0(alpha, gamma_nm)
alpha_nm = ballistic.alpha_nm0(alpha, gamma_nm)
ems_nm = ballistic.ems_nm0(ems, gamma_nm)

ene = np.arange(-0.1, 0.2, 0.005)
dens = np.empty_like(ene)
dens2 = np.empty_like(ene)
dens3 = np.empty_like(ene)

for i, e in enumerate(ene):
    dens[i] = ballistic.density1d_parabollic00(e, Enm, ems_nm, temperature)
    dens2[i] = ballistic.density1d_nonpara00(
        e, Enm, alpha_nm, ems_nm, temperature)
    dens3[i] = ballistic.density1d_all00(e, alpha, ems, temperature,
                                         W1, W2, 2, 2)

plt.plot(ene, dens, ene, dens2, ene, dens3)
plt.show()
