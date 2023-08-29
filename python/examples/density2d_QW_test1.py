#!/usr/bin/env python
# -*- coding: utf-8 -*-

# test of density2d_QW

import fetmodel
import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.constants as const
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
    ems = 0.067
    tOX = 20e-9
    temperature = 300
    W1 = args.width*1e-9
    W2 = 8e-9  # does not matter but set value for clarity
    alpha = fetmodel.alpha_NP(Eg, ems)
    # alpha_D = 0
    # alpha_G = 1
    alpha = fetmodel.alpha_NP(Eg, ems)
    Cox = epsOX*const.epsilon_0/tOX
    Cc = math.sqrt(const.elementary_charge*1e21*epsS*const.epsilon_0/(2*1))
    print(Cox, Cc)
    p = fetmodel.parameters_ballistic(alpha=alpha,
                                      Ceff=Cox*Cc/(Cox+Cc),
                                      ems=ems,
                                      W1=W1,
                                      W2=W2,
                                      nmax=args.nmax,
                                      mmax=1)
    p.output()

    EFermi_list = np.arange(-0.15, 0.1, 0.005)
    ns = np.empty_like(EFermi_list)
    ns0 = np.empty_like(EFermi_list)
    for i, EF0 in enumerate(EFermi_list):
        ns0[i] = fetmodel.density2d(EF0, p)
        ns[i] = fetmodel.density2d_QW_all(EF0, p)

    plt.plot(EFermi_list, ns, label='QW('+str(args.width)+'nm)')
    plt.plot(EFermi_list, ns0, label='no confinement (single subband)')
    plt.xlabel('Fermi Energy (eV)')
    plt.ylabel('2DEG density (m^-2)')
    plt.legend(loc='best')
    plt.show()
