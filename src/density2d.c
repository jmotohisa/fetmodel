/*
 *  density2d.c - Time-stamp: <Tue Aug 29 07:18:08 JST 2023>
 *
 *   Copyright (c) 2019  jmotohisa (Junichi Motohisa)  <motohisa@ist.hokudai.ac.jp>
 *
 *   Redistribution and use in source and binary forms, with or without
 *   modification, are permitted provided that the following conditions
 *   are met:
 *
 *   1. Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *
 *   2. Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 *   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 *   TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *  $Id: density2d.c 2019-08-02 21:18:09 jmotohisa $
 */

/*! 
  @file density2d.c 
  @brief 
  @author J. Motohisa
  @date
*/

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <tgmath.h>

#include <gsl/gsl_sf_fermi_dirac.h>
#include <gsl/gsl_const_mks.h>
#include <gsl/gsl_integration.h>

#include "ballistic_common.h"

#define GLOBAL_VALUE_DEFINE
#include "density2d.h"

/*!
  @brief
  @param[in]
  @param[out]
  @param[in,out]
  @return
*/

double density2d0(double EFermi,double Enm, double ems, double temp)
{
  double dos2D=MASS(ems)/(M_PI*GSL_CONST_MKS_PLANCKS_CONSTANT_HBAR*GSL_CONST_MKS_PLANCKS_CONSTANT_HBAR)*kBT0;
  return (dos2D*gsl_sf_fermi_dirac_0((EFermi-Enm)*BETA));
}

// parabolic band

// Quantization energy: rectangular Quantum Well
double Ep_n_rectQW(double ems, double W1, int n)
{
  double ene;
  ene=(GSL_CONST_MKS_PLANCKS_CONSTANT_HBAR*GSL_CONST_MKS_PLANCKS_CONSTANT_HBAR*M_PI*M_PI)/(2*MASS(ems)*GSL_CONST_MKS_ELECTRON_VOLT);
  return (ene*(n*n/(W1*W1)));
}

double density2d_QW0(double EFermi,double ems, double temp,
		     double W1, int n)
{
  double Enm = Ep_n_rectQW(ems, W1, n);
  return(density2d0(EFermi,Enm, ems, temp));
}

double density2d_QW_all0(double EFermi,double ems, double temp,
			 double W1, int nmax)
{
  int n;
  double sum;

  sum=0;
  for(n=1;n<=nmax;n++)
    sum += density2d_QW0(EFermi,ems,temp,W1,n);
  return(sum);
}

double density2d_QW_all(double EFermi, param_density2d_QW p)
{
  return density2d_QW_all0(EFermi,p.ems, p.temp,p.W1, p.nmax);
}

// nonparabolic band
// nonparabolicity parameter
// double Enm, double alpha_nm, double ems_nm
