/*
 *  pyfet.h - last saved: Time-stamp: <Tue Aug 20 10:17:07 JST 2019>
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
 *  $Id: pyfet.h 2019-07-12 14:25:09 jmotohisa $
 */

/*! 
  @file pyfet.h 
  @brief 
  @author J. Motohisa
*/

#ifndef _PYFET_H
#define _PYFET_H

#ifdef __cplusplus
extern "C" {
#endif

#ifdef	GLOBAL_VALUE_DEFINE
#define	GLOBAL
#else
#define	GLOBAL extern
#endif
  
  GLOBAL param_cMOSFET *param_cMOSFET_new(void );  
  GLOBAL param_cMESFET *param_cMESFET_new(void );  
  GLOBAL param_ballistic *param_ballistic_new(void );

  /* GLOBAL void Qapprox_cMOS_func(double *in_array,double *out_array, int size, param_cMOSFET p); */
  /* GLOBAL void Q_cMOS_func(double *in_array,double *out_array, int size, param_cMOSFET p); */

  /* GLOBAL void Ids0_cMOS_func(double *in_array,double *out_array, int size, double Vds,param_cMOSFET p); */
  /* GLOBAL void Ids_cMOS_func(double *in_array,double *out_array, int size, double Vds,param_cMOSFET p); */
  /* GLOBAL void Ids0_cMOS_R_func(double *in_array,double *out_array, int size, double Vds,param_cMOSFET p); */
  /* GLOBAL void Ids_cMOS_R_func(double *in_array,double *out_array, int size, double Vds,param_cMOSFET p); */
  
  /* GLOBAL void Ids_cMES_func(double *in_array,double *out_array, int size, double Vds,param_cMESFET p); */
  /* GLOBAL void Ids_cMES_R_func(double *in_array,double *out_array, int size, double Vds,param_cMESFET p); */

#undef GLOBAL_VALUE_DEFINE
#undef GLOBAL

#ifdef __cplusplus
}
#endif

#endif  // _PYFET_H
