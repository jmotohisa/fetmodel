/*
 *  mos1d.c - Time-stamp: <Sun Jun 22 19:09:06 JST 2025>
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
 *  $Id: mos1d.c 2020-01-08 20:08:38 jmotohisa $
 */

/*! 
  @file mos1d.c 
  @brief Change density, Drain current in 1D-MOS diode and MOSFET
  @author J. Motohisa
  @date
*/

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <tgmath.h>
#include <gsl/gsl_const_mks.h>

#include "config.h"
#include "ballistic_common.h"

#define GLOBAL_VALUE_DEFINE
#include "mos1d.h"

/*!
  @brief
  @param[in]
  @param[out]
  @param[in,out]
  @return
*/

/*!
  @brief
  @input betapsi beta * psi (beta = q/(kB*T), psi is surface potentiaol)
  @imput np0pp0 
  @return 
 */
double func_F(double betapsi,double np0pp0)
{
  return(sqrt(exp(-betapsi)+betapsi-1+np0pp0*(exp(betapsi)-betapsi-1)));
}

/*!
  @brief Space charge density in the inversion/accumulation layer
  @input psis surface potential
  @input p parameter of plMOSFET
  @return charge density
 */
double func_Qs_plMOS(double psis,param_plMOSFET p)
{
  double np0,pp0;
  np0=p.ni*p.ni/p.NA;
  pp0=p.NA;
  return func_QsMOS1D(psis,np0,pp0,p.eps_semi,p.temp);
}

/*!
  @brief Space charge density in the inversion/accumulation layer
  @input psis surface potential
  @input np0 minority carrier concentration
  @input pp0 majority carrier concentration
  @input eps_s dielectric constant of semiconductor
  @input temp temperature
  @return charge density
 */
double func_QsMOS1D(double psis,double np0,double pp0,double eps_s,double temp)
{
  double LD;
  LD=sqrt(kBT0*EPSILON(eps_s)/pp0)/GSL_CONST_MKS_ELECTRON_VOLT;
  return(sqrt(2.)*kBT0*EPSILON(eps_s)/(GSL_CONST_MKS_ELECTRON_VOLT*LD)*func_F(BETA*psis,np0/pp0));
}
