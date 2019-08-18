/*
 *  ballistic.h - last saved: Time-stamp: <Sat Aug 10 12:05:54 JST 2019>
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
 *  $Id: ballistic.h 2019-07-29 11:46:52 jmotohisa $
 */

/*! 
  @file ballistic.h 
  @brief 
  @author J. Motohisa
*/

#ifndef _BALLISTIC_H
#define _BALLISTIC_H

#ifdef __cplusplus
extern "C" {
#endif

#ifdef	GLOBAL_VALUE_DEFINE
#define	GLOBAL
#else
#define	GLOBAL extern
#endif

typedef struct param_ballistic_struct
{
  double EFermi;
  double VDS;
  double VGS;
  double alpha_D;
  double alpha_G;
  double Ceff;
  double alpha;
  double ems;
  double temp;
  double W1;
  double W2;
  int nmax;
  int mmax;
} param_ballistic;

typedef struct param_E0_struct
{
  double EFermi;
  double VDS;
  double VGS;
  double alpha_D;
  double alpha_G;
  double Ceff;
  param_density1d_all p_density1d_all;
} param_E0;
  
  GLOBAL double func_for_rootfind_E0_rect1d0(double ene0,double EFermi,
											 double VDS, double VGS,
											 double alpha_D, double alpha_G,
											 double Ceff,
											 double alpha, double ems, double temp,
											 double W1, double W2, int nmax, int mmax);
  GLOBAL double func_for_rootfind_E0_rect1d(double ene0,param_E0 *p);

  GLOBAL double E0_rect1d_root0(double EFermi,
					   double VDS, double VGS, double alpha_D, double alpha_G, double Ceff,
					   double alpha, double ems, double temp,
								double W1, double W2, int nmax, int mmax);
  GLOBAL double E0_rect1d_root(param_ballistic p);

 
#undef GLOBAL_VALUE_DEFINE
#undef GLOBAL

#ifdef __cplusplus
}
#endif

#endif  // _BALLISTIC_H
