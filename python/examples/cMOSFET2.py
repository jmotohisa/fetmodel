#!/usr/bin/env python
# -*- coding: utf-8 -*-

import fetmodel
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import scipy.constants as const
from scipy import integrate
import cMOSFET

if __name__ == '__main__':
    p = fetmodel.param_cMOSFET()
    p.radius = 50e-9
    p.Lg = 1e-6
    p.eps_semi = 15.15
    p.Rs = 0
    p.Rd = 0
    p.temp = 300
    p.ni = 1e21
    p.dphi = 0
    p.tox = 3e-9
    p.eps_ox = 20
    p.mue = 0.1
    p.Cox = p.eps_ox*8.85e-12/(p.radius*math.log(1+p.tox/p.radius))
    # p.Cox = fetmodel.Cox_radial(p.eps_ox,p.txo,p.radius)

    qq = np.arange(-20, -1,0.01)
    qfunc1=np.empty_like(qq)
    qfunc2=np.empty_like(qq)
    V=0
    Vgs=0.25
    
    # print(qfunc_cMOSFET(qq,V,Vgs,p))
    # print(qroot_brent(V,Vgs,p))
    # print(fetmodel.func_rootfind_Q_cMOSFET(qq,V,Vgs,p))
    
    # for i, q0 in enumerate(qq):
    #    qfunc1[i] = cMOSFET.qfunc_cMOSFET(10.**q0,V,Vgs,p)
    #    qfunc2[i] = fetmodel.func_rootfind_Q_cMOSFET(10.**q0,V,Vgs,p)

    # print(qfunc1-qfunc2)
    # plt.plot(qq,qfunc1,qq,qfunc2)
    # plt.show()

    ps=fetmodel.param_solver()

    Vds=0.5
    Ids_cutoff=1e-4
    Voff=cMOSFET.Voff_from_Ioff_cMOSFET(p, Vds, Ids_cutoff, ps)
    print(Voff)
    print(fetmodel.func_Ids_cMOSFET(Vds,Voff+0.25,p,ps))

    print(fetmodel.func_Qcharge2_cMOSFET(Vds, Voff, p))
    print(fetmodel.func_Qcharge2_cMOSFET(Vds, Voff+0.5, p))

    for i, q0 in enumerate(qq):
       qfunc1[i] = cMOSFET.qfunc_cMOSFET(10.**q0,Vds,Voff,p)
       qfunc2[i] = fetmodel.func_rootfind_Q_cMOSFET(10.**q0,Vds,Voff+0.25,p)
    
    # print(qfunc1-qfunc2)
    plt.plot(qq,qfunc1,qq,qfunc2)
    plt.show()
    
    print(cMOSFET.qroot_brent(Vds,Voff,p))
    print(cMOSFET.qroot_brent(Vds,Voff+0.25,p))
