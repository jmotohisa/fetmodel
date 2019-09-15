#!/usr/bin/env python
# -*- coding: utf-8 -*-

import fetmodel
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import scipy.constants as const
import ballistic1d


if __name__ == '__main__':
    p = fetmodel.param_ballistic_new()
    p.W1=5e-9
    p.W2=8e-9
    p.Eg=3.4
    p.ems=1/(1+20/p.Eg)
    p.alpha=(1-p.ems)**2/p.Eg
    n=2
    m=3
    Enmp_f = fetmodel.Ep_nm_rect1d(p.ems, p.W1, p.W2, n, m)
    gamma_nm_f = fetmodel.gamma_nm_NP(Enmp_f, p.alpha)
    Enm_f = fetmodel.E_nm_NP(p.alpha, gamma_nm_f)

    Enmp_b = ballistic1d.Ep_nm_rect1d(p.ems, p.W1, p.W2, n, m)
    gamma_nm_b = ballistic1d.gamma_nm_NP(Enmp_b, p.alpha)
    Enm_b = ballistic1d.E_nm_NP(p.alpha, gamma_nm_f)

    print(Enmp_f,Enmp_b)
    print(gamma_nm_f,gamma_nm_b)
    print(Enm_f,Enm_b)
