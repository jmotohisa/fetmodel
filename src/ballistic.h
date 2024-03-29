/*
 *  ballistic.h - last saved: Time-stamp: <Tue Aug 29 19:47:26 JST 2023>
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

  GLOBAL double func_for_findroot_E0_rect1dNP0(double ene0,double EFermi,
					       double VDS, double VGS,
					       double alpha_D, double alpha_G,
					       double Ceff,
					       double alpha, double ems, double temp,
					       double W1, double W2, int nmax, int mmax);
  
  GLOBAL double E0_rect1dNP_root0(double EFermi,
				  double VDS, double VGS,
				  double alpha_D, double alpha_G, double Ceff,
				  double alpha, double ems, double temp,
				  double W1, double W2, int nmax, int mmax);
  
  GLOBAL double Ids_ballistic1d_rect1dNP0(double VDS, double VGS, double EFs, 
					  double alpha_D, double alpha_G,
					  double Ceff,
					  double alpha, double ems, double temp,
					  double W1, double W2, int nmax, int mmax);
  
  // 2d (no confinement, single subband)
  GLOBAL double func_for_findroot_E0_2d0(double ene0,double EFermi,
					 double VDS, double VGS, 
					 double alpha_D, double alpha_G,
					 double Ceff,
					 double ems, double temp);
  GLOBAL double E0_2d_root0(double EFermi,
			    double VDS, double VGS, double alpha_D, double alpha_G, double Ceff,
			    double ems, double temp);
  GLOBAL double Ids_ballistic2d0_E0(double E0,
				    double VDS, double VGS, double EFs,
				    double alpha_D, double alpha_G,
				    double Ceff,
				    double ems, double temp);
  GLOBAL double Ids_ballistic2d0(double VDS, double VGS, double EFs,
				 double alpha_D, double alpha_G,
				 double Ceff,
				 double ems, double temp);

  // multiple subbands
  GLOBAL double func_for_findroot_E0_QW0(double ene0,double EFermi,
					 double VDS, double VGS,
					 double alpha_D, double alpha_G,
					 double Ceff,
					 double alpha, double ems, double temp,
					 double W1, int nmax);
  
  GLOBAL double E0_QW_root0(double EFermi,
			    double VDS, double VGS,
			    double alpha_D, double alpha_G, double Ceff,
			    double alpha, double ems, double temp,
			    double W1, int nmax);
  
  GLOBAL double Ids_ballistic2d_QW0(double VDS, double VGS, double EFs, 
				    double alpha_D, double alpha_G,
				    double Ceff,
				    double alpha, double ems, double temp,
				    double W1, int nmax);

#undef GLOBAL_VALUE_DEFINE
#undef GLOBAL

#ifdef __cplusplus
}
#endif

#endif  // _BALLISTIC_H
