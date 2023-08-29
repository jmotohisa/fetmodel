#!/usr/bin/env python
# -*- coding: utf-8 -*-

# test of ballisti 2d FET with QW confinement

import fetmodel
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import ballistic2d
import argparse


def get_args():
    parser = argparse.ArgumentParser(
        description='test of density2d')
    parser.add_argument('-w', '--width',
                        nargs='?',
                        type=float,
                        help='well width (nm)',
                        default=10)
    parser.add_argument('-n', '--nmax',
                        nargs='?',
                        type=int,
                        help='number of subbands',
                        default=10)
    args = parser.parse_args()

    return(args)


if __name__ == '__main__':
    args = get_args()

    Eg = 0.36
    epsOX = 8.5
    epsS = 8.9
    tOX = 20e-9
    temperature = 300
    ems = 0.067
    W1 = args.width*1e-9
    W2 = 8e-9  # does not matter but set value for clarity
    alpha = fetmodel.alpha_NP(Eg, ems)
    Cox = epsOX*8.85e-12/tOX
    Cc = math.sqrt(1.6e-19*1e21/(2*epsS*8.86e-12*1))
    alpha_D = 0
    alpha_G = 1

    p = fetmodel.param_ballistic()
    p.ems = ems
    p.alpha = alpha
    p.W1 = W1
    p.W2 = W2

    p.EFermi = -0.1
    p.VDS = 0
    p.VGS = 0
    p.alpha_D = alpha_D
    p.alpha_G = alpha_G
    p.Ceff = Cox*Cc/(Cox+Cc)
    p.temp = temperature
    p.nmax = args.nmax
    p.mmax = 2  # does not matter but set value for clarity

    Vgs = np.arange(-0.1, 1, 0.01)
    Ids1 = np.empty_like(Vgs)
    Ids2 = np.empty_like(Vgs)
    Vds = 0.5

    for i, Vgs0 in enumerate(Vgs):
        Ids1[i] = fetmodel.Ids_ballistic2d_QW(Vds, Vgs0, p, 0)*1e-3
        Ids2[i] = ballistic2d.Ids_ballistic2d(Vds, Vgs0, p, 0)*1e-3

    fig, ax = plt.subplots()
    ax.plot(Vgs, Ids1, label='QW (width='+str(args.width)+' nm)')
    ax.plot(Vgs, Ids2, label='single subband')
    ax.set_yscale("log")
    plt.xlabel('Gate Voltage (V)')
    plt.ylabel('Drain Current (mA/um)')
    plt.legend(loc='best')
    plt.show()
