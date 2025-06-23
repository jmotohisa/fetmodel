#!/usr/bin/env python
# -*- coding: utf-8 -*-

import scipy.constants as const
import fetmodel


class parameters_ballistic:

    def __init__(self,
                 alpha_D=0.,
                 alpha_G=1.,
                 Ceff=-1.,
                 alpha=0.,
                 ems=0.2,
                 temp=300.,
                 W1=5e-9,
                 W2=8e-9,
                 nmax=1,
                 mmax=1):

        self.alpha_D = alpha_D
        self.alpha_G = alpha_G
        self.Ceff = Ceff
        self.alpha = alpha
        self.ems = ems
        self.temp = temp
        self.W1 = W1
        self.W2 = W2
        self.nmax = nmax
        self.mmax = mmax

    def output(self):
        print("alpha_D=", self.alpha_D)
        print("alpha_G=", self.alpha_G)
        print("Ceff=", self.Ceff)
        print("alpha=", self.alpha)
        print("ems=", self.ems)
        print("temp=", self.temp)
        print("W1=", self.W1)
        print("W2=", self.W2)
        print("nmax=", self.nmax)
        print("mmax=", self.mmax)


class parameters_cMOSFET:
    def __init__(self,
                 radius=6.25e-9,
                 Lg=1.5e-6,
                 eps_semi=11.6,
                 Rs=0.,
                 Rd=0.,
                 Cox=-1.,
                 temp=300.,
                 ni=1.45e16,
                 dphi=0.,
                 tox=1.5e-9,
                 eps_ox=3.9,
                 mue=0.04):

        self.cMOSFET = fetmodel.param_cMOSFET()
        self.cMOSFET.radius = radius
        self.cMOSFET.Lg = Lg
        self.cMOSFET.eps_semi = eps_semi
        self.cMOSFET.Rs = Rs
        self.cMOSFET.Rd = Rd
        if Cox < 0:
            self.cMOSFET.Cox = fetmodel.Cox_radial_area(eps_ox, tox, radius)
        else:
            self.cMOSFET.Cox = Cox
        self.cMOSFET.temp = temp
        self.cMOSFET.ni = ni
        self.cMOSFET.dphi = dphi
        self.cMOSFET.tox = tox
        self.cMOSFET.eps_ox = eps_ox
        self.cMOSFET.mue = mue

    def output(self):
        print("radius=", self.cMOSFET.radius)
        print("Lg=", self.cMOSFET.Lg)
        print("eps_semi=", self.cMOSFET.eps_semi)
        print("Rs=", self.cMOSFET.Rs)
        print("Rd=", self.cMOSFET.Rd)
        print("Cox=", self.cMOSFET.Cox)
        print("temp=", self.cMOSFET.temp)
        print("ni=", self.cMOSFET.ni)
        print("dphi=", self.cMOSFET.dphi)
        print("tox=", self.cMOSFET.tox)
        print("eps_ox=", self.cMOSFET.eps_ox)
        print("mue=", self.cMOSFET.mue)


class parameters_cMESFET:
    def __init__(self,
                 radius=6.25e-9,
                 Lg=1.5e-6,
                 eps_semi=11.6,
                 Rs=0.,
                 Rd=0.,
                 temp=300.,
                 mue=0.04,
                 Nd=1e22,
                 Vbi=0.7,):

        self.cMESFET = fetmodel.param_cMESFET()
        self.cMESFET.radius = radius
        self.cMESFET.Lg = Lg
        self.cMESFET.eps_semi = eps_semi
        self.cMESFET.Rs = Rs
        self.cMESFET.Rd = Rd
        self.cMESFET.temp = temp
        self.cMESFET.mue = mue
        self.cMESFET.Nd = Nd
        self.cMESFET.Vbi = Vbi

    def output(self):
        print("radius=", self.cMESFET.radius)
        print("Lg=", self.cMESFET.Lg)
        print("eps_semi=", self.cMESFET.eps_semi)
        print("Rs=", self.cMESFET.Rs)
        print("Rd=", self.cMESFET.Rd)
        print("temp=", self.cMESFET.temp)
        print("mue=", self.cMESFET.mue)
        print("Nd=", self.cMESFET.Nd)
        print("Vbi=", self.cMESFET.Vbi)


class parameters_plMOSFET:
    def __init__(self,
                 Lg=1.5e-6,
                 eps_semi=11.6,
                 Rs=0.,
                 Rd=0.,
                 Cox=-1.,
                 temp=300.,
                 ni=1.45e16,
                 dphi=0.,
                 tox=1.5e-9,
                 eps_ox=3.9,
                 mue=0.04,
                 NA=1e21):

        self.plMOSFET = fetmodel.param_plMOSFET()
        self.plMOSFET.Lg = Lg
        self.plMOSFET.eps_semi = eps_semi
        self.plMOSFET.Rs = Rs
        self.plMOSFET.Rd = Rd
        if Cox < 0:
            self.plMOSFET.Cox = const.epsilon_0*eps_ox/tox
        else:
            self.plMOSFET.Cox = Cox
        self.plMOSFET.temp = temp
        self.plMOSFET.ni = ni
        self.plMOSFET.dphi = dphi
        self.plMOSFET.tox = tox
        self.plMOSFET.eps_ox = eps_ox
        self.plMOSFET.mue = mue
        self.plMOSFET.NA = NA

    def output(self):
        print("Lg=", self.plMOSFET.Lg)
        print("eps_semi=", self.plMOSFET.eps_semi)
        print("Rs=", self.plMOSFET.Rs)
        print("Rd=", self.plMOSFET.Rd)
        print("Cox=", self.plMOSFET.Cox)
        print("temp=", self.plMOSFET.temp)
        print("ni=", self.plMOSFET.ni)
        print("dphi=", self.plMOSFET.dphi)
        print("tox=", self.plMOSFET.tox)
        print("eps_ox=", self.plMOSFET.eps_ox)
        print("mue=", self.plMOSFET.mue)
        print("NA=", self.plMOSFET.NA)


class parameters_solver:

    def __init__(self,
                 left=0.,
                 right=1.
                 ):

        self.param_solver = fetmodel.param_solver()
        self.param_solver.left = left
        self.param_solver.right = right

    def output(self):
        print("left=", self.param_solver.left)
        print("right=", self.param_solver.right)
