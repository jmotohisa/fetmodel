/*
 *  cfet.c - Time-stamp: <Fri Jul 12 15:34:54 JST 2019>
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
 *  $Id: cfet.c 2019-07-12 14:25:01 jmotohisa $
 */

/*! 
  @file cfet.c 
  @brief 
  @author J. Motohisa
  @date
*/

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <tgmath.h>
#include "../ctl-io.h"

#define GLOBAL_VALUE_DEFINE
#include "cfet.h"

/*!
  @brief
  @param[in]
  @param[out]
  @param[in,out]
  @return
*/

param_cMOSFET *param_cMOSFET_new()
{
  param_cMOSFET *p = malloc(sizeof(param_cMOSFET));
  if(p ==NULL) {
	perror("malloc");
	exit(1);
  }
  return(p);
}

param_cMESFET *param_cMESFET_new()
{
  param_cMESFET *p = malloc(sizeof(param_cMESFET));
  if(p ==NULL) {
	perror("malloc");
	exit(1);
  }
  return(p);
}

void set_global_cMOSFET(param_cMOSFET p)
{
  radius=p.radius;
  Lg=p.Lg;
  eps_semi=p.eps_semi;
  Rs=p.Rs;
  Rd=p.Rd;
  Cox=p.Cox;
  temp=p.temp;
  ni=p.ni;
  dphi=p.dphi;
  tox=p.tox;
  eps_ox=p.eps_ox;
  mue=p.mue;
}

double Ids_cMOS(double Vds, double Vgs,param_cMOSFET p)
{
  set_global_cMOSFET(p);
  return(Ids_cMOSFET(Vds,Vgs));
}
	
  
