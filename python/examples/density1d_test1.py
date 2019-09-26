#!/usr/bin/env python
# -*- coding: utf-8 -*-

# test of density1d

import fetmodel
import numpy as np
import matplotlib.pyplot as plt

Eg = 3.4
ems = 0.2
alpha = fetmodel.alpha_NP(Eg, ems)
temperature = 300

W1 = 5e-9
W2 = 8e-9

Enmp = fetmodel.Ep_nm_rect1d(ems, W1, W2, 1, 1)
gamma_nm = fetmodel.gamma_nm_NP(Enmp, alpha)
Enm = fetmodel.E_nm_NP(alpha, gamma_nm)
alpha_nm = fetmodel.alpha_nm_NP(alpha, gamma_nm)
ems_nm = fetmodel.ems_nm_NP(ems, gamma_nm)

print(alpha, gamma_nm, Enm, alpha_nm, ems_nm)

ene = np.arange(-0.1, 0.2, 0.005)
dens1 = np.empty_like(ene)
dens2 = np.empty_like(ene)
dens3 = np.empty_like(ene)

for i, e in enumerate(ene):
    dens1[i] = fetmodel.density1d0(e, Enm, ems_nm, temperature)
    dens2[i] = fetmodel.density1d_NP0(
        e, Enm, alpha_nm, ems_nm, temperature)
    dens3[i] = fetmodel.density1d_rect1dNP_all0(e, alpha, ems, temperature,
                                                W1, W2, 2, 3)

plt.plot(ene, dens1,label="parabolic")
plt.plot(ene, dens2,label="nonpalabolic")
plt.plot(ene, dens3,label="multiband")
plt.legend(loc='best')
plt.show()

dens1 = np.empty_like(ene)
dens2 = np.empty_like(ene)
dens3 = np.empty_like(ene)
for i, e in enumerate(ene):
    dens1[i] = fetmodel.density1d_rect1dNP_all0(e, alpha, ems, temperature,
                                                5e-9, 8e-9, 2, 3)
    dens2[i] = fetmodel.density1d_rect1dNP_all0(e, alpha, ems, temperature,
                                                5e-9, 16e-9, 2, 3)
    dens3[i] = fetmodel.density1d_rect1dNP_all0(e, alpha, ems, temperature,
                                                5e-9, 24e-9, 2, 3)
    
plt.plot(ene, dens1,label="W2=8nm")
plt.plot(ene, dens2,label="W2=16nm")
plt.plot(ene, dens3,label="W2=24nm")
plt.legend(loc='best')
plt.xlim([-0.1,0.2])
plt.ylim([0,1.5e9])
plt.show()
