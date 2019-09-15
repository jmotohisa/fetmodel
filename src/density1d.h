/*
 *  density1d.h - last saved: Time-stamp: <Thu Sep 12 14:13:34 JST 2019>
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
 *  $Id: density1d.h 2019-07-29 11:32:49 jmotohisa $
 */

/*! 
  @file density1d.h 
  @brief 
  @author J. Motohisa
*/

#ifndef _DENSITY1D_H
#define _DENSITY1D_H

#ifdef __cplusplus
extern "C" {
#endif

#ifdef	GLOBAL_VALUE_DEFINE
#define	GLOBAL
#else
#define	GLOBAL extern
#endif

  typedef struct param_density1d_struct
  {
	double temp;
	double EFermi;
	double Enm;
	double alphanm;
	double emsnm;
  } param_density1d;
  
  typedef struct param_density1d_rect_struct
  {
	double alpha;
	double ems;
	double temp;
	double W1;
	double W2;
	int nmax;
	int mmax;
  } param_density1d_rect;

  // parabolic band  
  GLOBAL double density1d0(double EFermi, double Enm, double ems,double temp);
  GLOBAL double density1d_rect1d0(double EFermi, double ems,double temp,
								  double W1, double W2, int n, int m);
  GLOBAL double density1d_rect1d_all0(double EFermi, double ems, double temp,
									  double W1, double W2, int nmax, int mmax);
  GLOBAL double density1d_rect1d_all(double EFermi,param_density1d_rect p);

  // nonparabolic band
  GLOBAL double density1d_NP0(double EFermi, double Enm, double alpha_nm, double ems_nm, double temp);
  GLOBAL double density1d_rect1dNP0(double EFermi,double alphaNP, double ems, double temp,
									double W1, double W2, int n, int m);
  GLOBAL double density1d_rect1dNP_all0(double EFermi,double alphaNP, double ems, double temp,
										double W1, double W2, int nmax, int mmax);

  GLOBAL double density1d_NP(param_density1d params);
  GLOBAL double density1d_rect1dNP_all(double EFermi,param_density1d_rect p);

  // utility functions
  GLOBAL double Ep_nm_rect1d(double ems, double W1, double W2, int n , int m);
  GLOBAL double Ep_n_radial1d(double ems,double radius,int n);
  GLOBAL double Ep_nm_radial1d(double ems,double radius,int n, int m);

  GLOBAL double alpha_NP(double Eg, double ems);
  GLOBAL double E_nm_NP(double alphaNP, double gamma_nm);
  GLOBAL double alpha_nm_NP(double alphaNP, double gamma_nm);
  GLOBAL double ems_nm_NP(double ems, double gamma_nm);
  GLOBAL double gamma_nm_NP(double Enm, double alphaNP);

  GLOBAL double gamma_nm_rect1dNP(double alphaNP,double ems,double W1,double W2,int n,int m);
  GLOBAL double E_nm_rect1dNP(double alphaNP,double ems,double W1,double W2,int n,int m);
  
#undef GLOBAL_VALUE_DEFINE
#undef GLOBAL

#ifdef __cplusplus
}
#endif

#endif  // _DENSITY1D_H
