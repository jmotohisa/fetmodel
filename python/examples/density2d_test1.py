#!/usr/bin/env python
# -*- coding: utf-8 -*-

# test of density2d

import fetmodel
import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.constants as const

if __name__ == '__main__':
    Eg = 0.36
    epsOX = 8.5
    epsS = 8.9
    ems = 0.067
    tOX = 20e-9
    temperature = 300
    W1 = 10e-9
    W2 = 8e-9
    alpha = fetmodel.alpha_NP(Eg, ems)
    # alpha_D = 0
    # alpha_G = 1
    alpha = fetmodel.alpha_NP(Eg, ems)
    Cox = epsOX*const.epsilon_0/tOX
    Cc = math.sqrt(const.elementary_charge*1e21*epsS*const.epsilon_0/(2*1))
    print(Cox,Cc)
    p=fetmodel.parameters_ballistic(alpha=alpha,
                                    Ceff=Cox*Cc/(Cox+Cc),
                                    ems=ems,
                                    W1=W1,
                                    W2=W2,
                                    nmax=1,
                                    mmax=1)
    p.output()

    EFermi_list = np.arange(-0.15, 0.1, 0.005)
    ns = np.empty_like(EFermi_list)
    for i, EF0 in enumerate(EFermi_list):
        ns[i] = fetmodel.density2d0(EF0, 0, p.ems, p.temp)
        
    plt.plot(EFermi_list, ns)
    plt.xlabel('Fermi Energy (eV)')
    plt.ylabel('2DEG density (m^-2)')
    plt.show()
        
