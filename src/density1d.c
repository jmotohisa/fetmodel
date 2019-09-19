/*
 *  density1d.c - Time-stamp: <Mon Sep 16 21:55:49 JST 2019>
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
 *  $Id: density1d.c 2019-07-29 11:32:18 jmotohisa $
 */

/*! 
  @file density1d.c 
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
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>

#include "ballistic_common.h"

#define GLOBAL_VALUE_DEFINE
#include "density1d.h"

#include "fd.h"

double density1d_NP(param_density1d params);
double density1d_rect1d_all(double EFermi,param_density1d_rect p);
double density1d_rect1dNP_all(double EFermi,param_density1d_rect p);

/*!
  @brief
  @param[in]
  @param[out]
  @param[in,out]
  @return
*/

// parabollic band

// Quantization energy: rectangular QWR
double Ep_nm_rect1d(double ems, double W1, double W2, int n , int m)
{
  double ene;
  ene=(GSL_CONST_MKS_PLANCKS_CONSTANT_HBAR*GSL_CONST_MKS_PLANCKS_CONSTANT_HBAR*M_PI*M_PI)/(2*MASS(ems)*GSL_CONST_MKS_ELECTRON_VOLT);
  return(ene*(n*n/(W1*W1)+m*m/(W2*W2)));
}

double Ep_n_radial1d(double ems,double radius,int n)
{
  double ene;
  double jzero=gsl_sf_bessel_zero_J0(n);
  ene=(GSL_CONST_MKS_PLANCKS_CONSTANT_HBAR*GSL_CONST_MKS_PLANCKS_CONSTANT_HBAR)/(2*MASS(ems)*GSL_CONST_MKS_ELECTRON_VOLT);
  return(ene*(jzero*jzero)/(radius*radius));
}

double Ep_nm_radial1d(double ems,double radius,int n, int m)
{
  double ene;
  double jzero=gsl_sf_bessel_zero_Jnu ((double) m, n);
  ene=(GSL_CONST_MKS_PLANCKS_CONSTANT_HBAR*GSL_CONST_MKS_PLANCKS_CONSTANT_HBAR)/(2*MASS(ems)*GSL_CONST_MKS_ELECTRON_VOLT);
  return(ene*(jzero*jzero)/(radius*radius));
}

double density1d0(double EFermi, double Enm, double ems,double temp)
{
  double d0=sqrt(2*MASS(ems)*kBT0/M_PI)/GSL_CONST_MKS_PLANCKS_CONSTANT_HBAR;
  double ene=(EFermi-Enm)/kBT;
  double retval;
  retval=fd_mhalf(ene);
  /* retval=gsl_sf_fermi_dirac_mhalf (ene); */
  return (d0*retval);
}

double density1d_rect1d0(double EFermi, double ems,double temp,
						 double W1, double W2, int n, int m)
{
  double Enm = Ep_nm_rect1d(ems, W1, W2, n , m);
  return(density1d0(EFermi, Enm, ems,temp));
}

double density1d_rect1d_all0(double EFermi, double ems, double temp,
							   double W1, double W2, int nmax, int mmax)
{
  int n,m;
  double sum;

  sum=0;
  for(n=1;n<=nmax;n++)
	for(m=1;m<=mmax;m++)
	  {
		sum += density1d_rect1d0(EFermi,ems,temp,W1,W2,n,m);
	  }
  return(sum);
}

double density1d_rect1d_all(double EFermi,param_density1d_rect p)
{
  return(density1d_rect1d_all0(EFermi,p.ems,p.temp,p.W1,p.W2,p.nmax,p.mmax));
}

// nonparabolic band
// nonparabolicity parameter
// double Enm, double alpha_nm, double ems_nm

double alpha_NP(double Eg, double ems)
{
  return((1-ems)*(1-ems)/Eg);
}

double E_nm_NP(double alphaNP, double gamma_nm)
{
  return((gamma_nm-1)/(2*alphaNP));
}

double alpha_nm_NP(double alphaNP, double gamma_nm)
{
  return(alphaNP/gamma_nm);
}

double ems_nm_NP(double ems, double gamma_nm)
{
  return(ems*gamma_nm);
}

double gamma_nm_NP(double Enm, double alphaNP)
{
  return(sqrt(1+4*alphaNP*Enm));
}

// rectangular QWR, nonparaboic band
// quantization energy with nonparaboic band correction
double gamma_nm_rect1dNP(double alphaNP,double ems,double W1,double W2,int n,int m)
{
  double Ep_nm=Ep_nm_rect1d(ems,W1,W2, n,m);
  return(sqrt(1+4*alphaNP*Ep_nm));
}											

double E_nm_rect1dNP(double alphaNP,double ems,double W1,double W2,int n,int m)
{
  double gamma_nm=gamma_nm_rect1dNP(alphaNP,ems,W1,W2,n,m);
  return(E_nm_NP(alphaNP,gamma_nm));
}

double alpha_nm_rect1dNP(double alpha,double ems,double W1,double W2,int n,int m)
{
  double gamma_nm=gamma_nm_rect1dNP(alpha,ems,W1,W2,n,m);
  return(alpha_nm_NP(alpha,gamma_nm));
}

double ems_nm_rect1dNP(double alpha,double ems,double W1,double W2,int n,int m)
{
  double gamma_nm=gamma_nm_rect1dNP(alpha,ems,W1,W2,n,m);
  return(ems_nm_NP(ems,gamma_nm));
}

// electron density
// reduced density of states:
double dos1D0_red(double ene,double alpha)
{
  if(ene>0)
	return((1+2*alpha*ene)/sqrt(ene*(1+alpha*ene)));
  else
	return 0;
}

double func_for_integ_density_1d0(double eta,double EFermi,double Enm,double emsnm,double alphanm,double temp)
{
  double eta0=(EFermi-Enm)/(kBT);
  return (dos1D0_red(eta,alphanm*kBT)/(1+exp(eta-eta0)));
}

double func_for_integ_density_1d(double ene, void *params)
{
  param_density1d *p=(param_density1d *) params;
  return(func_for_integ_density_1d0(ene,p->EFermi,p->Enm, p->emsnm, p->alphanm,p->temp));
}

// test of integration
/* double f (double x, void * params) { */
/*   double alpha = *(double *) params; */
/*   double f = log(alpha*x) / sqrt(x); */
/*   return f; */
/* } */

/* int gsl_test() */
/* { */
/*   gsl_integration_workspace * w = */
/* 	gsl_integration_workspace_alloc(1000); */
/*   double result, error; */
/*   double expected = -4.0; */
/*   double alpha = 1.0; */
/*   gsl_function F; */
/*   F.function = &f; */
/*   F.params = &alpha; */
/*   gsl_integration_qags(&F, 0, 1, 0, 1e-7, 1000, */
/* 					   w, &result, &error); */
/* } */

double density1d_NP(param_density1d params)
{
  gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
  gsl_function F;
  double epsabs=1e-7;
  double epsrel=1e-7;
  size_t limit=1000;
  double result, abserr;
  double temp=params.temp;
  int status;
  gsl_error_handler_t old_handler,qagiu_handler;
  
  F.function=&func_for_integ_density_1d;
  F.params = &params;
  /* int gsl_integration_qagiu (gsl_function * f, double a, double epsabs, double epsrel, size_t limit, gsl_integration_workspace * workspace, double *result, double *abserr); */
  
  gsl_set_error_handler_off();
  /* old_handler = gsl_set_error_handler(&qagiu_handler); */
  
  status=gsl_integration_qagiu(&F,0,epsabs,epsrel,limit,w,&result,&abserr);
  gsl_integration_workspace_free(w);
  /* printf("%d\n",status); */
  if(status!=GSL_SUCCESS) {
	printf("Error in density1d_NP:status=%d\n",status);
	GSL_ERROR_VAL("argument lies on singularity",
				  GSL_ERANGE, GSL_NAN);
	}
  gsl_set_error_handler(NULL);
  double d0=sqrt(2*MASS(params.emsnm)*kBT0)/(GSL_CONST_MKS_PLANCKS_CONSTANT_HBAR*M_PI);
  /* gsl_set_error_handler(old_handler); */
  return(d0*result);
}
				  
double density1d_NP0(double EFermi, double Enm, double alpha_nm, double ems_nm, double temp)
{
  param_density1d p;
  p.temp=temp;
  p.EFermi=EFermi;
  p.Enm=Enm;
  p.alphanm=alpha_nm;
  p.emsnm=ems_nm;
  
  return (density1d_NP(p));
}

double density1d_rect1dNP0(double EFermi,double alphaNP, double ems, double temp,
						   double W1, double W2, int n, int m)
{
  double gamma_nm, Enm; //,alpha_nm, ems_nm;

  /* alphaNP = alpha_NP(Eg, ems); */
  gamma_nm = gamma_nm_rect1dNP(alphaNP, ems, W1, W2, n, m);
  Enm = E_nm_NP(alphaNP, gamma_nm);
  return(density1d_NP0(EFermi, Enm,
					   alpha_nm_NP(alphaNP, gamma_nm),
					   ems_nm_NP(ems, gamma_nm),
					   temp));
}

double density1d_rect1dNP_all0(double EFermi,double alphaNP, double ems, double temp,
							   double W1, double W2, int nmax, int mmax)
{
  int n,m;
  double sum;

  /* alphaNP = alpha_NP(Eg, ems); */
  sum=0;
  for(n=1;n<=nmax;n++)
	for(m=1;m<=mmax;m++)
	  {
		sum += density1d_rect1dNP0(EFermi,alphaNP, ems, temp,W1,W2,n,m);
	  }
  return(sum);
}

double density1d_rect1dNP_all(double EFermi,param_density1d_rect p)
{
  return
	(density1d_rect1dNP_all0(EFermi,p.alpha, p.ems, p.temp, p.W1, p.W2, p.nmax, p.mmax));
}
  
