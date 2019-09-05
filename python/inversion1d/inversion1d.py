#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
inversion charge in classical MOS diode
"""

import numpy as np
import math
import scipy.constants as const
import matplotlib.pyplot as plt
import fetmodel


def func_Qcharge1dMOS0(psis, ni, eps_semi, Na, temp):
    beta = const.elementary_charge/(temp*const.Boltzmann)
    beta1 = 1/beta
    LD = math.sqrt(beta1*eps_semi*const.epsilon_0/(const.elementary_charge*Na))
    q0 = math.sqrt(2)*eps_semi*const.epsilon_0*temp * \
        const.Boltzmann/(const.elementary_charge*LD)
    pp0 = Na
    np0 = ni**2/Na
    f = (math.exp(-beta*psis)+beta*psis-1)+np0 / \
        pp0*(math.exp(beta*psis) - beta*psis-1)
    return math.sqrt(abs(f))*q0


if __name__ == '__main__':
    p = fetmodel.param_cMOSFET_new()
    p.eps_semi = 11.6
    p.temp = 300.
    p.ni = 1.45e16
    Na = 1e21

    psis = np.arange(-0.5, 1.2, 0.002)
    qcharge1 = np.empty_like(psis)
    qcharge2 = np.empty_like(psis)
    for i, psis0 in enumerate(psis):
        qcharge1[i] = func_Qcharge1dMOS0(
            psis0, p.ni, p.eps_semi, 1e21, p.temp)*1e-4
        qcharge2[i] = func_Qcharge1dMOS0(
            psis0, p.ni, p.eps_semi, 1e20, p.temp)*1e-4

    fig, ax = plt.subplots()
    ax.plot(psis, qcharge1, label='Na=1e15cm^-3')
    ax.plot(psis, qcharge2, label='Na=1e14cm^-3')
    ax.set_yscale("log")
    ax.set_ylim([1e-9, 1e-4])
    plt.ylabel('Qs (C/cm^2)')
    plt.xlabel('psi_s (eV)')
    plt.legend(loc='best')
    plt.show()
