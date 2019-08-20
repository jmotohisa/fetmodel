/*
 *  ballistic.c - Time-stamp: <Wed Aug 21 07:03:47 JST 2019>
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
#include <gsl/gsl_sf_fermi_dirac.h>

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
  param_density1d_rect p_density1d_rect;
  
  p_density1d_rect.alpha=alpha;
  p_density1d_rect.ems=ems;
  p_density1d_rect.temp=temp;
  p_density1d_rect.W1=W1;
  p_density1d_rect.W2=W2;
  p_density1d_rect.nmax=nmax;
  p_density1d_rect.mmax=mmax;

  n1d_S=density1d_rect1d_all(EFermi-ene0,    p_density1d_rect);
  n1d_D=density1d_rect1d_all(EFermi-ene0-VDS,p_density1d_rect);
  return(alpha_D*VDS + alpha_G*VGS - (n1d_S + n1d_D)/(2*Ceff)*GSL_CONST_MKS_ELECTRON_VOLT);
}

double func_for_findroot_E0_rect1d(double ene0,param_E0 *p)
{
  return(func_for_findroot_E0_rect1d0(ene0,
									  p->EFermi, p->VDS, p->VGS,
									  p->alpha_D, p->alpha_G,
									  p->Ceff,
									  p->p_density1d_rect.alpha, p->p_density1d_rect.ems,
									  p->p_density1d_rect.temp,
									  p->p_density1d_rect.W1, p->p_density1d_rect.W2,
									  p->p_density1d_rect.nmax, p->p_density1d_rect.mmax));
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
  param_density1d_rect p_density1d_rect;
  param_E0 params;
  
  p_density1d_rect.alpha=alpha;
  p_density1d_rect.ems=ems;
  p_density1d_rect.temp=temp;
  p_density1d_rect.W1=W1;
  p_density1d_rect.W2=W2;
  p_density1d_rect.nmax=nmax;
  p_density1d_rect.mmax=mmax;
  
  params.p_density1d_rect = p_density1d_rect;
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

// current

double Ids_ballistic1d_recdt1dNP0(double VDS, double VGS,
								double EFs, double EFermi,
								double alpha_D, double alpha_G,
								double Ceff,
								double alpha, double ems, double temp,
								double W1, double W2, int nmax, int mmax)
{
  double Enm,E0;
  double ids1,ids2;
  int n,m;
  ids1=0;
  ids2=0;
  E0=E0_rect1d_root0(EFermi,VDS, VGS, alpha_D, alpha_G, Ceff,
					 alpha, ems, temp,
					 W1, W2, nmax, nmax);
  
  for(n=1;n<=nmax;n++)
    for(m=1;m<=mmax;m++)
	  {
		Enm=Ep_nm_rect1d(ems,W1,W2,n,m);
		ids1 += gsl_sf_fermi_dirac_0(BETA*(EFs-Enm-E0));
		ids2 += gsl_sf_fermi_dirac_0(BETA*(EFs-Enm-E0-VDS));
	  }
  
  return(2*(ids1-ids2)*GSL_CONST_MKS_ELECTRON_VOLT/GSL_CONST_MKS_PLANCKS_CONSTANT_H*kBT0);
}

double Ids_ballistic1d_recdt1dNP(param_ballistic p,double EFs)
{
  return(Ids_ballistic1d_recdt1dNP0(p.VDS,p.VGS,EFs,p.EFermi,
									p.alpha_D, p.alpha_G, p.Ceff,
									p.alpha, p.ems, p.temp,
									p.W1, p.W2, p.nmax, p.nmax));
}
