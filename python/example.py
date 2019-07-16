import pycfet
import math
# import numpy as np
# import matplotlib.pyplot as plt
# import h5py

p = pycfet.param_cMOSFET_new()
p.radius = 50e-9
p.Lg = 1e-6
p.eps_semi = 11.6
p.Rs = 0
p.Rd = 0
p.Cox = 0.0233
p.temp = 300
p.ni = 1.45e16
p.dphi = 0
p.tox = 1.5e-9
p.eps_ox = 3.9
p.mue = 0.04

Vgs = 1
Vds = 1

print(pycfet.Qapprox_cMOS(Vgs, p))
print(pycfet.Qapprox_cMOS(Vgs-Vds, p))
print(pycfet.Q_cMOS(Vgs, p))
print(pycfet.Q_cMOS(Vgs-Vds, p))
print(pycfet.Ids_cMOS(Vds, Vgs, p))