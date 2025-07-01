#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import math
import scipy.constants as const
import matplotlib.pyplot as plt
import fetmodel

if __name__ == '__main__':
    pSi = fetmodel.parameters_plMOSFET()
    pSi.eps_semi = 11.6
    pSi.temp = 300.
    pSi.ni = 1.45e16
    pSi.NA = 1e22
    pInAs = fetmodel.parameters_plMOSFET()
    pInAs.eps_semi = 15.15
    pInAs.temp = 300.
    pInAs.ni = 1e21
    pInAs.NA = 1e22

    psis = np.arange(-0.5, 1.2, 0.002)
    qcharge1 = np.empty_like(psis)
    qcharge2 = np.empty_like(psis)
    for i, psis0 in enumerate(psis):
        # qcharge1[i] = fetmodel.func_Qcharge1dMOS0(
        #     psis0, pSi.ni, pSi.eps_semi, 1e22, pSi.temp)*1e-4
        # qcharge2[i] = fetmodel.func_Qcharge1dMOS0(
        #     psis0, pInAs.ni, pInAs.eps_semi, 1e22, pInAs.temp)*1e-4
        qcharge1[i] = fetmodel.Q_plMOS(psis0, pSi)*1e-4
        qcharge2[i] = fetmodel.Q_plMOS(psis0, pInAs)*1e-4

    # qcharge1 = (qcharge2/qcharge1)
    fig, ax = plt.subplots()
    ax.plot(psis, qcharge1, label='Si')
    ax.plot(psis, qcharge2, label='InAs')
    ax.set_yscale("log")
    ax.set_ylim([1e-9, 1e-4])
    plt.ylabel(r'Qs (C/cm^2)')
    plt.xlabel(r'$\psi_s$ (eV)')
    plt.legend(loc='best')
    plt.show()
