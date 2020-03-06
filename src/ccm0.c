/*
 *  ccm0.c - Time-stamp: <Thu Mar 05 21:44:59 JST 2020>
 *
 *   Copyright (c) 2020  jmotohisa (Junichi Motohisa)  <motohisa@ist.hokudai.ac.jp>
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
 *  $Id: ccm0.c 2020-03-05 21:35:21 jmotohisa $
 */

/*! 
  @file ccm0.c 
  @brief 
  @author J. Motohisa
  @date
*/

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <tgmath.h>
#include <float.h>
#include "ccm.h"
#include "PhysicalConstants.h"

#define GLOBAL_VALUE_DEFINE
#include "ccm0.h"

/*!
  @brief
  @param[in]
  @param[out]
  @param[in,out]
  @return
*/

#define POW2(X) (X)*(X)
#define kBT (KBC*p.temp)
#define Vth (kBT/EC)
#define Q0 (4*p.eps_semi*EPSILON/p.radius*Vth)
#define delta (EC*EC*p.ni/(kBT*p.eps_semi*EPSILON))
#define V0 (p.dphi + Vth*log(8/(delta*p.radius*p.radius)))

// approximate form of Q (eq.17 in Iniguez Trans. ED)
double Q_approx_0(double V, double Vgs,param_cMOSFET p)
{
  double qp0,Vt0,dVt0;
  qp0 = Q_approx0_0(V,Vgs,V0,0,p);
  Vt0=V0 + 2*Vth*log(1+qp0/Q0);
  dVt0 = (2*p.Cox*POW2(Vth)/Q0)*qp0/(Q0+qp0);
  return(Q_approx0_0(V,Vgs,Vt0,dVt0,p));
}

// unified approximated formula for Q (eq. 13 in Iniguez et al., Trans. ED)
double Q_approx0_0(double V,double Vgs,double Vt,double deltaVt,
				   param_cMOSFET p)
{
  double a,b,c1,c2;
  a=2*p.Cox*Vth*Vth/Q0;
  b=2*Vth*log(1+exp((Vgs-Vt+deltaVt-V)/(2*Vth)));
  c1=sqrt(a*a+b*b)-a;
  c2=b*b/(2*a);
  if(c1 < DBL_EPSILON)
	return(p.Cox*c2);
  else
	return(p.Cox*c1);
}
