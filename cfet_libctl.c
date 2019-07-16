/*
 *  cfet_libctl.c - Time-stamp: <Tue Jul 16 14:42:46 JST 2019>
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
 *  $Id: cfet_libctl.c 2019-07-12 16:07:51 jmotohisa $
 */

/*! 
  @file cfet_libctl.c 
  @brief 
  @author J. Motohisa
  @date
*/

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <tgmath.h>

#include "ctl-io.h"
#include "PhysicalConstants.h"
#include "ccm.h"

//#define GLOBAL_VALUE_DEFINE
//#include "cfet_libctl.h"


/*!
  @brief
  @param[in]
  @param[out]
  @param[in,out]
  @return
*/

// current

void get_global_cMOSFET(param_cMOSFET *p)
{
  p->radius=radius;
  p->Lg=Lg;
  p->eps_semi=eps_semi;
  p->Rs=Rs;
  p->Rd=Rd;
  p->Cox=Cox;
  p->temp=temp;
  p->ni=ni;
  p->dphi=dphi;
  p->tox=tox;
  p->eps_ox=eps_ox;
  p->mue=mue;
}

void get_global_cMESFET(param_cMESFET *p)
{
  p->radius=radius;
  p->Lg=Lg;
  p->eps_semi=eps_semi;
  p->Rs=Rs;
  p->Rd=Rd;
  /* p->Cox=Cox; */
  /* p->temp=temp; */
  /* p->ni=ni; */
  /* p->dphi=dphi; */
  /* p->tox=tox; */
  /* p->eps_ox=eps_ox; */
  /* p->mue=mue; */
  p->Nd=Nd;
  p->Vbi=Vbi;
}

number Ids_cMOSFET(number Vds,number Vgs)
{
  param_cMOSFET p;
  get_global_cMOSFET(&p);
  return(Ids0_cMOSFET(Vds,Vgs,p));
}

number Qapprox_cMOSFET(number Vgs)
{
  param_cMOSFET p;
  get_global_cMOSFET(&p);
  return(Qapprox_cMOS0(Vgs,p));
}

number Q_cMOSFET(number Vgs)
{
  param_cMOSFET p;
  get_global_cMOSFET(&p);
  return(Q_cMOS0(Vgs,p));
}

number Ids_cMOSFET_R(number Vds,number Vgs)
{
  param_cMOSFET p;
  get_global_cMOSFET(&p);
  return(Ids0_cMOSFET_R(Vds,Vgs,p));
}

number Ids_cMESFET(number Vds,number Vgs)
{
  param_cMESFET p;
  get_global_cMESFET(&p);
  return(Ids0_cMESFET(Vds,Vgs,p));
}

number Ids_cMESFET_R(number Vds,number Vgs)
{
  param_cMESFET p;
  get_global_cMESFET(&p);
  return(Ids0_cMESFET_R(Vds,Vgs,p));
}
