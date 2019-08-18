/*
 *  ballistic.c - Time-stamp: <Sat Aug 10 19:43:21 JST 2019>
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
 *  $Id: ballistic.c 2019-07-29 09:31:40 jmotohisa $
 */

/*! 
  @file ballistic.c 
  @brief 
  @author J. Motohisa
  @date
*/

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <tgmath.h>
#include <gsl/gsl_const_mks.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#include "ballistic_common.h"
#include "density1d.h"

#define GLOBAL_VALUE_DEFINE
#include "ballistic.h"


/*!
  @brief
  @param[in]
  @param[out]
  @param[in,out]
  @return
*/

// 1D, rectangular cross section, with nonparabolicity

double func_for_findroot_E0_rect1d0(double ene0,double EFermi,
									double VDS, double VGS,
									double alpha_D, double alpha_G,
									double Ceff,
									double alpha, double ems, double temp,
									double W1, double W2, int nmax, int mmax)					
{
  double n1d_S,n1d_D;
  param_density1d_all p_density1d_all;
  
  p_density1d_all.alpha=alpha;
  p_density1d_all.ems=ems;
  p_density1d_all.temp=temp;
  p_density1d_all.W1=W1;
  p_density1d_all.W2=W2;
  p_density1d_all.nmax=nmax;
  p_density1d_all.mmax=mmax;

  n1d_S=density1d_all0(EFermi-ene0,    p_density1d_all);
  n1d_D=density1d_all0(EFermi-ene0-VDS,p_density1d_all);
  return(alpha_D*VDS + alpha_D*VGS - (n1d_S + n1d_D)/(2*Ceff)*GSL_CONST_MKS_ELECTRON_VOLT);
}

double func_for_findroot_E0_rect1d(double ene0,param_E0 *p)
{
  return(func_for_findroot_E0_rect1d0(ene0,
									  p->EFermi, p->VDS, p->VGS,
									  p->alpha_D, p->alpha_G,
									  p->Ceff,
									  p->p_density1d_all.alpha, p->p_density1d_all.ems,
									  p->p_density1d_all.temp,
									  p->p_density1d_all.W1, p->p_density1d_all.W2,
									  p->p_density1d_all.nmax, p->p_density1d_all.mmax));
  /* double n1d_S,n1d_D; */
  
  /* n1d_S=density1d_all0(p->EFermi-ene0,       p->p_density1d_all); */
  /* n1d_D=density1d_all0(p->EFermi-ene0-p->VDS,p->p_density1d_all); */
  /* return(p->alpha_D*p->VDS + p->alpha_D*p->VGS - (n1d_S + n1d_D)/(2*p->Ceff)*GSL_CONST_MKS_ELECTRON_VOLT); */
}

double func_for_findroot_E0_rect1d00(double ene0,void *pp)
{
  param_E0 *p = (param_E0 *) pp;
  double e0=func_for_findroot_E0_rect1d(ene0,p);
  return(ene0+e0);
}

double E0_rect1d_root_brent(param_E0 params,
							double low, double high)
{
  gsl_function F;
  int status;
  int iter = 0, max_iter = 100;
  double r;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;

  F.function = &func_for_findroot_E0_rect1d00;
  F.params = &params;

  //FindRoots/Q/L=(low) qfunc_cMOSFET,param_cMOSFET;

  T= gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc(T);

  gsl_root_fsolver_set(s,&F,low,high);
  do {
        iter++;
        status = gsl_root_fsolver_iterate(s);
        r=gsl_root_fsolver_root(s);
        low = gsl_root_fsolver_x_lower(s);
        high = gsl_root_fsolver_x_upper(s);
        status = gsl_root_test_interval(low,high,0,0.001);
  } while (status==GSL_CONTINUE && iter < max_iter);

  return(r);
}

double E0_rect1d_root0(double EFermi,
					   double VDS, double VGS, double alpha_D, double alpha_G, double Ceff,
					   double alpha, double ems, double temp,
					   double W1, double W2, int nmax, int mmax)
{
  double low,high;
  param_density1d_all p_density1d_all;
  param_E0 params;
  
  p_density1d_all.alpha=alpha;
  p_density1d_all.ems=ems;
  p_density1d_all.temp=temp;
  p_density1d_all.W1=W1;
  p_density1d_all.W2=W2;
  p_density1d_all.nmax=nmax;
  p_density1d_all.mmax=mmax;
  
  params.p_density1d_all = p_density1d_all;
  params.EFermi = EFermi;
  params.VDS = VDS;
  params.VGS = VGS;
  params.alpha_D = alpha_D;
  params.alpha_G = alpha_G;
  params.Ceff = Ceff;
  low=-1;
  high=1;
  return(E0_rect1d_root_brent(params, low, high));
  
}

double E0_rect1d_root(param_ballistic p)
{
  return(E0_rect1d_root0(p.EFermi,p.VDS, p.VGS, p.alpha_D, p.alpha_G, p.Ceff,
						 p.alpha, p.ems, p.temp,
						 p.W1, p.W2, p.nmax, p.nmax));
}
