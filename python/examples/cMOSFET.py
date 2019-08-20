#!/usr/bin/env python
# -*- coding: utf-8 -*-

# sample Python script of cyrlindircal MOSFET
# B. Iniguez et al., IEEE Trans. Elec. Dev. 52, No.8, Auguust 2015 (p.1868)

import pyfet
import math
import numpy as np
import matplotlib.pyplot as plt
# import h5py

p = pyfet.param_cMOSFET_new()
p.radius = 6.25e-9
p.Lg = 1e-6
p.eps_semi = 11.6
p.Rs = 0
p.Rd = 0
p.temp = 300
p.ni = 1.45e16
p.dphi = 0
p.tox = 1.5e-9
p.eps_ox = 3.9
p.mue = 0.04
p.Cox = p.eps_ox*8.85e-12/(p.radius*math.log(1+p.tox/p.radius))
# p.Cox = pyfet.Cox_radial(p.eps_ox,p.txo,p.radius)

Vgs = np.arange(0, 2, 0.01)

# Fig.2
Q1 = np.empty_like(Vgs)
Q2 = np.empty_like(Vgs)
pyfet.Qapprox_cMOS(Vgs, Q1, p)
pyfet.Q_cMOS(Vgs, Q2, p)
# Q1 = pyfet.Qapprox_cMOS(Vgs, p)
# Q2 = pyfet.Q_cMOS(Vgs, p)

plt.plot(Vgs, Q1)
plt.plot(Vgs, Q2)
plt.show()

# Fig.3
Ids01 = np.empty_like(Vgs)
Ids02 = np.empty_like(Vgs)
pyfet.Ids_cMOS(Vgs, Ids01, 0.1, p)
pyfet.Ids0_cMOS(Vgs, Ids02, 0.1, p)
plt.plot(Vgs, Ids01)
plt.plot(Vgs, Ids02)
plt.show()

# Fig.4
pyfet.Ids_cMOS(Vgs, Ids01, 1, p)
pyfet.Ids0_cMOS(Vgs, Ids02, 1, p)
plt.plot(Vgs, Ids01)
plt.plot(Vgs, Ids02)
plt.show()
