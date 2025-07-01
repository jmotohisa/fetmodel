#!/usr/bin/env python
# -*- coding: utf-8 -*-

import fetmodel
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import scipy.constants as const
from scipy import integrate
# import h5py

p = fetmodel.parameters_plMOSFET()

p.radius = 6.25e-9
p.Lg = 1e-6
p.eps_semi = 11.6
p.Rs = 0
p.Rd = 0
p.temp = 300
p.ni = 1.45e16
p.dphi = 0
p.tox = 20e-9
p.eps_ox = 3.9
p.mue = 0.04
p.Cox = p.eps_ox*const.epsilon_0/p.tox  # planar mosfet
p.NA = 1e22

# p.Cox = fetmodel.Cox_radial(p.eps_ox,p.txo,p.radius)


# psiB=psiB_func(p,NA)

NA = 1e22
Vbs = np.linspace(-10, 0)
Vth1 = np.empty_like(Vbs)
for i, Vbs0 in enumerate(Vbs):
    Vth1[i] = fetmodel.Vth_body(Vbs0, p)

p.NA = 3e21
Vth2 = np.empty_like(Vbs)
for i, Vbs0 in enumerate(Vbs):
    Vth2[i] = fetmodel.Vth_body(Vbs0, p)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlim([0, 10])
ax.set_xlabel(r'$-V_{BS}$ [V]', size=18)
ax.set_ylim([0.7, 1.8])
ax.set_ylabel(r'Threshold Voltage [V]', size=18)
ax.plot(-Vbs, Vth1, label=r'$N_A=1\times 10^{16} \mathrm{cm}^{-3}$')
ax.plot(-Vbs, Vth2, label=r'$N_A=3\times 10^{15} \mathrm{cm}^{-3}$')
ax.legend(loc='best', fontsize=18)
ax.tick_params(labelsize=18)
plt.show()

Vds = np.linspace(0, 5)
Vgs = np.linspace(1, 5, 5)

Ids1 = np.zeros([Vgs.shape[0], Vds.shape[0]])
Ids2 = np.zeros([Vgs.shape[0], Vds.shape[0]])
Ids3 = np.zeros([Vgs.shape[0], Vds.shape[0]])

p.output()
print(Vgs)

for i, Vgs0 in enumerate(Vgs):
    for j, Vds0 in enumerate(Vds):
        Ids1[i][j] = fetmodel.Ids1_plMOSFET(Vgs0, Vds0, p)
        Ids2[i][j] = fetmodel.Ids0_plMOSFET(
            Vgs0, Vds0, fetmodel.Vth_plMOSFET(p), p)
        Ids3[i][j] = fetmodel.Ids2_plMOSFET(Vgs0, Vds0, p)

fig = plt.figure()
ax = fig.add_subplot(111)
for j, ids0 in enumerate(Ids2):
    if (j == 0):
        ax.plot(Vds, Ids2[j, :], color='black', label=r'Eq.(3.2)')
        ax.plot(Vds, Ids1[j, :], color='red', label=r'Eq.(3.3)')
        ax.plot(Vds, Ids3[j, :], color='blue', label=r'Eq.(3.4)')
    else:
        ax.plot(Vds, Ids2[j, :], color='black')
        ax.plot(Vds, Ids1[j, :], color='red')
        ax.plot(Vds, Ids3[j, :], color='blue')

ax.set_xlim([0, 5])
ax.set_xlabel(r'$V_{DS}$ [V]', size=18)
# ax.set_ylim([0,2])
ax.set_ylabel(r'$I_{DS}$ [A/m]', size=18)
# ax.plot(Vds, Vth2,label=r'$N_A=3\times 10^{15} \mathrm{cm}^{-3}$')
ax.legend(loc='best', fontsize=18)
ax.tick_params(labelsize=18)
plt.show()

Vgs = np.linspace(1, 5)
Vds = 5
Ids1 = np.empty_like(Vgs)
Ids2 = np.empty_like(Vgs)
Ids3 = np.empty_like(Vgs)

for i, Vgs0 in enumerate(Vgs):
    Ids1[i] = fetmodel.Ids1_plMOSFET(Vgs0, Vds, p)
    Ids2[i] = fetmodel.Ids0_plMOSFET(
        Vgs0, Vds, fetmodel.Vth_plMOSFET(p), p)
    Ids3[i] = fetmodel.Ids2_plMOSFET(Vgs0, Vds, p)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlim([0, 5])
ax.set_xlabel(r'$V_{GS}$ [V]', size=18)
# ax.set_ylim([0,2])
ax.set_ylabel(r'$I_{DS}$ [A/m]', size=18)
ax.plot(Vgs, Ids1, label=r'$N_A=3\times 10^{15} \mathrm{cm}^{-3}$')
ax.plot(Vgs, Ids2, label=r'charge sheet model')
ax.plot(Vgs, Ids3, label=r'second order')
# ax.plot(Vds, Vth2,label=r'$N_A=3\times 10^{15} \mathrm{cm}^{-3}$')
ax.legend(loc='best', fontsize=18)
ax.tick_params(labelsize=18)
plt.show()

psiS = np.linspace(0, 3)
res = np.empty_like(psiS)
for i, psiS0 in enumerate(psiS):
    res[i] = fetmodel.func_psiS_findroot1_plMOSFET(psiS0, p, 1)

p.NA = 1e23
p.tox = 10e-9
p.Cox = const.epsilon_0*p.eps_ox/p.tox

Vgs = 1
Vds = np.linspace(0, 3)
psiS1 = np.empty_like(Vds)
psiS2 = np.empty_like(Vds)
psiS3 = np.empty_like(Vds)
psiS4 = np.empty_like(Vds)
for i, Vds0 in enumerate(Vds):
    Vgs = 1
    psiS1[i] = fetmodel.find_psiS_plMOSFET(p, Vgs, Vds0)
    Vgs = 3
    psiS2[i] = fetmodel.find_psiS_plMOSFET(p, Vgs, Vds0)
    Vgs = 5
    psiS3[i] = fetmodel.find_psiS_plMOSFET(p, Vgs, Vds0)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(Vds, psiS1)
ax.plot(Vds, psiS2)
ax.plot(Vds, psiS3)
psiS4 = 2*fetmodel.psiB_func(p)+Vds
ax.plot(Vds, psiS4, '.')
plt.show()

Vds = np.linspace(0, 5)
Vgs = np.linspace(1, 5, 5)
print(Vgs)

Ids1 = np.zeros([Vgs.shape[0], Vds.shape[0]])
Ids2 = np.zeros([Vgs.shape[0], Vds.shape[0]])
Ids3 = np.zeros([Vgs.shape[0], Vds.shape[0]])
Ids4 = np.zeros([Vgs.shape[0], Vds.shape[0]])

for i, Vgs0 in enumerate(Vgs):
    for j, Vds0 in enumerate(Vds):
        Ids1[i][j] = fetmodel.Ids1_plMOSFET(Vgs0, Vds0, p)
        Ids2[i][j] = fetmodel.Ids0_plMOSFET(
            Vgs0, Vds0, fetmodel.Vth_plMOSFET(p), p)
        Ids3[i][j] = fetmodel.Ids2_plMOSFET(Vgs0, Vds0, p)
        Ids4[i][j] = fetmodel.Ids_plMOSFET(Vgs0, Vds0, p)

fig = plt.figure()
ax = fig.add_subplot(111)
for j, ids0 in enumerate(Ids2):
    if (j == 0):
        ax.plot(Vds, Ids2[j, :], color='black', label=r'Eq.(3.2)')
        ax.plot(Vds, Ids1[j, :], color='red', label=r'Eq.(3.3)')
        ax.plot(Vds, Ids3[j, :], color='blue', label=r'Eq.(3.4)')
        ax.plot(Vds, Ids4[j, :], color='orange', label=r'strict')
    else:
        ax.plot(Vds, Ids2[j, :], color='black')
        ax.plot(Vds, Ids1[j, :], color='red')
        ax.plot(Vds, Ids3[j, :], color='blue')
        ax.plot(Vds, Ids4[j, :], color='orange')

ax.set_xlim([0, 5])
ax.set_xlabel(r'$V_{DS}$ [V]', size=18)
# ax.set_ylim([0,2])
ax.set_ylabel(r'$I_{DS}$ [A/m]', size=18)
# ax.plot(Vds, Vth2,label=r'$N_A=3\times 10^{15} \mathrm{cm}^{-3}$')
ax.legend(loc='best', fontsize=18)
ax.tick_params(labelsize=18)
