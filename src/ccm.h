/*
 *  ccm.h - last saved: Time-stamp: <Mon Mar 09 15:31:31 JST 2020>
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
	/* double Eg; */
  } param_cMOSFET;
  
  typedef struct
  {
	double	radius;
	double	Lg;
	double	eps_semi;
	double	Rs;
	double	Rd;
	/* double	Cox; */
	double	temp;
	/* double	ni; */
	/* double	dphi; */
	/* double	tox; */
	/* double	eps_ox; */
	double	mue;
  	double	Nd;
	double	Vbi;
  } param_cMESFET;

  typedef struct
  {
	double left;
	double right;
  } param_solver;

  GLOBAL double func_rootfind_Q_cMOSFET(double qq,double V, double Vgs, param_cMOSFET p);

  GLOBAL double func_Qcharge_cMOSFET(double V, double Vgs, param_cMOSFET p, param_solver ps);
  GLOBAL double func_Qcharge2_cMOSFET(double V, double Vgs, param_cMOSFET p);
  GLOBAL double func_rootfind_Q_cMOSFET(double qq, double V, double Vgs, param_cMOSFET p);
  GLOBAL double func_Ids_cMOSFET(double Vds,double Vgs,param_cMOSFET p,param_solver ps);
  GLOBAL double func_Ids2_cMOSFET(double Vds,double Vgs,param_cMOSFET p);
  GLOBAL double func_Ids_cMOSFET_R(double Vds,double Vgs,param_cMOSFET p,param_solver ps);
  GLOBAL double func_Ids2_cMOSFET_R(double Vds,double Vgs,param_cMOSFET p);

  GLOBAL double func_Ids_cMESFET(double Vds,double Vgs,param_cMESFET p, param_solver ps);
  GLOBAL double func_Ids_cMESFET_R(double Vds,double Vgs,param_cMESFET p, param_solver ps);

#undef GLOBAL_VALUE_DEFINE
#undef GLOBAL

#ifdef __cplusplus
}
#endif

#endif  // _CCM_H
