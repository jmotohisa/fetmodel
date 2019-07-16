/*
 *  ccm.h - last saved: Time-stamp: <Tue Jul 16 14:36:18 JST 2019>
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
 *  $Id: ccm.h 2019-07-12 16:06:01 jmotohisa $
 */

/*! 
  @file ccm.h 
  @brief 
  @author J. Motohisa
*/

#ifndef _CCM_H
#define _CCM_H

#ifdef __cplusplus
extern "C" {
#endif

#ifdef	GLOBAL_VALUE_DEFINE
#define	GLOBAL
#else
#define	GLOBAL extern
#endif

  typedef struct
  {
	double	radius;
	double	Lg;
	double	eps_semi;
	double	Rs;
	double	Rd;
	double	Cox;
	double	temp;
	double	ni;
	double	dphi;
	double	tox;
	double	eps_ox;
	double	mue;
  } param_cMOSFET;
  
  typedef struct
  {
	double	radius;
	double	Lg;
	double	eps_semi;
	double	Rs;
	double	Rd;
	/* double	Cox; */
	/* double	temp; */
	/* double	ni; */
	/* double	dphi; */
	/* double	tox; */
	/* double	eps_ox; */
	/* double	mue; */
  	double	Nd;
	double	Vbi;
  } param_cMESFET;

  GLOBAL double Ids0_cMOSFET(double Vds,double Vgs,param_cMOSFET p);
  GLOBAL double Ids0_cMESFET(double Vds,double Vgs,param_cMESFET p);
  GLOBAL double Ids0_cMOSFET_R(double Vds,double Vgs,param_cMOSFET cMOS);
  GLOBAL double Ids0_cMESFET_R(double Vds,double Vgs,param_cMESFET cMES);

  GLOBAL double Qapprox_cMOS0(double Vgs, param_cMOSFET p);
  GLOBAL double Q_cMOS0(double Vgs, param_cMOSFET p);

#undef GLOBAL_VALUE_DEFINE
#undef GLOBAL

#ifdef __cplusplus
}
#endif

#endif  // _CCM_H
