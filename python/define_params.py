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

        self.alpha_D=alpha_D
        self.alpha_G=alpha_G
        self.Ceff=Ceff
        self.alpha=alpha
        self.ems=ems
        self.temp=temp
        self.W1=W1
        self.W2=W2
        self.nmax=nmax
        self.mmax=nmax

    def output(self):
        print("alpha_D=",self.alpha_D)
        print("alpha_G=",self.alpha_G)
        print("Ceff=",self.Ceff)
        print("alpha=",self.alpha)
        print("ems=",self.ems)
        print("temp=",self.temp)
        print("W1=",self.W1)
        print("W2=",self.W2)
        print("nmax=",self.nmax)
        print("mmax=",self.nmax)

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

        self.radius=radius
        self.Lg=Lg
        self.eps_semi=eps_semi
        self.Rs=Rs
        self.Rd=Rd
        if Cox<0 :
            self.Cox=fetmodel.Cox_radial_area(eps_ox,tox,radius)
        else:
            self.Cox=Cox
        self.temp=temp
        self.ni=ni
        self.dphi=dphi
        self.tox=tox
        self.eps_ox=eps_ox
        self.mue=mue
        
    def output(self):
        print("Lg=",self.Lg)
        print("eps_semi=",self.eps_semi)
        print("Rs=",self.Rs)
        print("Rd=",self.Rd)
        print("Cox=",self.Cox)
        print("temp=",self.temp)
        print("ni=",self.ni)
        print("dphi=",self.dphi)
        print("tox=",self.tox)
        print("eps_ox=",self.eps_ox)
        print("mue=",self.mue)
