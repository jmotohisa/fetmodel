/*
 *  density1d.c - Time-stamp: <Sat Aug 10 12:05:06 JST 2019>
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

#include "ballistic_common.h"

#define GLOBAL_VALUE_DEFINE
#include "density1d.h"

/*!
  @brief
  @param[in]
  @param[out]
  @param[in,out]
  @return
*/

double alphaNP00(double Eg, double ems)
{
  return((1-ems)*(1-ems)/Eg);
}

double Ep_nm00(double ems, double W1, double W2, int n , int m)
{
  double ene;
  ene=(GSL_CONST_MKS_PLANCKS_CONSTANT_HBAR*GSL_CONST_MKS_PLANCKS_CONSTANT_HBAR*M_PI*M_PI)/(2*MASS(ems)*GSL_CONST_MKS_ELECTRON_VOLT)
	;
  return(ene*(n*n/(W1*W1)+m*m/(W2*W2)));
}

double gamma_nm00(double alpha,double ems,double W1,double W2,int n,int m)
{
  double Ep_nm=Ep_nm00(ems,W1,W2, n,m);
  return(sqrt(1+4*alpha*Ep_nm));
}											

// nonparabolicity parameter
// double Enm, double alpha_nm, double ems_nm

double E_nm0(double alpha, double gamma_nm)
{
  return((gamma_nm-1)/(2*alpha));
}

double alpha_nm0(double alpha, double gamma_nm)
{
  return(alpha/gamma_nm);
}

double ems_nm0(double ems, double gamma_nm)
{
  return(ems*gamma_nm);
}

double E_nm00(double alpha,double ems,double W1,double W2,int n,int m)
{
  double gamma_nm=gamma_nm00(alpha,ems,W1,W2,n,m);
  return(E_nm0(alpha,gamma_nm));
}

double alpha_nm00(double alpha,double ems,double W1,double W2,int n,int m)
{
  double gamma_nm=gamma_nm00(alpha,ems,W1,W2,n,m);
  return(alpha_nm0(alpha,gamma_nm));
}

double ems_nm00(double alpha,double ems,double W1,double W2,int n,int m)
{
  double gamma_nm=gamma_nm00(alpha,ems,W1,W2,n,m);
  return(ems_nm0(ems,gamma_nm));
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

double func_density_1d0(double eta,double EFermi,double Enm,double emsnm,double alphanm,double temp)
{
  double eta0=(EFermi-Enm)/(kBT);
  return (dos1D0_red(eta,alphanm*kBT)/(1+exp(eta-eta0)));
}

double func_density_1d(double ene, void *params)
{
  param_density1d *p=(param_density1d *) params;
  return(func_density_1d0(ene,p->EFermi,p->Enm, p->emsnm, p->alphanm,p->temp));
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

double density1d_nonpara0(param_density1d params)
{
  gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
  gsl_function F;
  double epsabs=1e-7;
  double epsrel=1e-7;
  size_t limit=1000;
  double result, abserr;

  F.function=&func_density_1d;
  F.params = &params;
/* int gsl_integration_qagiu (gsl_function * f, double a, double epsabs, double epsrel, size_t limit, gsl_integration_workspace * workspace, double *result, double *abserr); */

  gsl_integration_qagiu(&F,0,epsabs,epsrel,limit,w,&result,&abserr);

  gsl_integration_workspace_free(w);
  return(result);
}
				  
double density1d_nonpara00(double EFermi, double Enm, double alpha_nm, double ems_nm, double temp)
{
  param_density1d p;
  p.temp=temp;
  p.EFermi=EFermi;
  p.Enm=Enm;
  p.alphanm=alpha_nm;
  p.emsnm=ems_nm;
  
  double d0=sqrt(2*MASS(p.emsnm)*kBT0)/(GSL_CONST_MKS_PLANCKS_CONSTANT_HBAR*M_PI);
  return (d0*density1d_nonpara0(p));
}

double density1d_parabollic00(double EFermi, double Enm, double ems,double temp)
{
  double d0=sqrt(2*MASS(ems)*kBT0/M_PI)/GSL_CONST_MKS_PLANCKS_CONSTANT_HBAR;
  return (d0*gsl_sf_fermi_dirac_mhalf ((EFermi-Enm)/kBT));
}

double density1d_all00(double EFermi,double alpha, double ems, double temp,
					   double W1, double W2, int nmax, int mmax)
{
  int n,m;
  double sum;
  double gamma_nm, Enm,alpha_nm, ems_nm;

  /* alpha = alphaNP00(Eg, ems); */
  sum=0;
  for(n=1;n<=nmax;n++)
	for(m=1;m<=mmax;m++)
	  {
		gamma_nm = gamma_nm00(alpha, ems, W1, W2, n, m);
		Enm = E_nm0(alpha, gamma_nm);
		alpha_nm = alpha_nm0(alpha, gamma_nm);
		ems_nm = ems_nm0(ems, gamma_nm);
		sum = sum+density1d_nonpara00(EFermi, Enm, alpha_nm, ems_nm, temp);
	  }
  return(sum);
}

double density1d_all0(double EFermi,param_density1d_all p)
{
  return
	(density1d_all00(EFermi,p.alpha, p.ems, p.temp, p.W1, p.W2, p.nmax, p.mmax));
}
  
