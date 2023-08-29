/*
 *  density2d.h - last saved: Time-stamp: <Tue Aug 29 07:16:35 JST 2023>
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
 *  $Id: density2d.h 2019-08-10 12:15:22 jmotohisa $
 */

/*! 
  @file density2d.h 
  @brief 
  @author J. Motohisa
*/

#ifndef _DENSITY2D_H
#define _DENSITY2D_H

#ifdef __cplusplus
extern "C" {
#endif

#ifdef	GLOBAL_VALUE_DEFINE
#define	GLOBAL
#else
#define	GLOBAL extern
#endif

    typedef struct param_density2d_struct
  {
	double temp;
	double EFermi;
	double Enm;
	double alphanm;
	double emsnm;
  } param_density2d;
  
  typedef struct param_density2d_QW_struct
  {
	double alpha;
	double ems;
	double temp;
	double W1;
	int nmax;
  } param_density2d_QW;


  GLOBAL double density2d0(double EFermi,double Enm, double ems, double temp);
  GLOBAL double density2d_QW0(double EFermi,double ems, double temp,double W1, int n);
  GLOBAL double density2d_QW_all0(double EFermi,double ems, double temp,
				  double W1, int nmax);

#undef GLOBAL_VALUE_DEFINE
#undef GLOBAL

#ifdef __cplusplus
}
#endif

#endif  // _DENSITY2D_H
