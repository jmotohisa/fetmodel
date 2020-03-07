/*
 *  cfet_libctl.c - Time-stamp: <Sun Mar 08 06:52:50 JST 2020>
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
#include "../src/ccm.h"
#include "../src/capacitor.h"

#define GLOBAL_VALUE_DEFINE
#include "cfet_libctl.h"

/*!
  @brief
  @param[in]
  @param[out]
  @param[in,out]
  @return
*/

double get_radius(params_NWFET p)
{
  double radius;
  switch(p.sizes.which_subclass) {
  case NW_RECT:
	radius = p.sizes.subclass.NW_rect_data->W1;
	break;
  case NW_RADIAL:
  default:
	radius = p.sizes.subclass.NW_radial_data->radius;
	break;
  }
  return radius;
}

double get_radius2(params_NWFET p)
{
  double radius;
  switch(p.sizes.which_subclass) {
  case NW_RECT:
	radius = p.sizes.subclass.NW_rect_data->W2;
	break;
  case NW_RADIAL:
  default:
	radius = p.sizes.subclass.NW_radial_data->radius;
	break;
  }
  return radius;
}

double get_Cox(params_NWFET p)
{
  double Cox0;
  params_cMOSFET p0;
  function *CoxFunc;
  SCM s;
  double x0=0;
  switch(p.which_subclass) {
  case PARAMS_CMOSFET:
  default:
	p0=*(p.subclass.params_cMOSFET_data);
	if(p0.use_Cox_funcp) {
	  CoxFunc = &(p0.Cox_func);
	  s=gh_call1(*CoxFunc,ctl_convert_number_to_scm(x0));
	  Cox0 = ctl_convert_number_to_c(s);
	}
	else
	  {
		switch(p.sizes.which_subclass) {
		case NW_RECT:
		  Cox0=Cox_rect_area(p.eps_ox,p.tox,get_radius(p),get_radius2(p));
		  break;
		case NW_RADIAL:
		default:
		  Cox0=Cox_radial_area(p.eps_ox,p.tox,get_radius(p));
		  //Cox_radial(double epsOX, double tOX, double radius);
		  break;
		}
	  }
	break;
  }
  return(Cox0);
}

// Get global and store local valuable 
void get_global_cMOSFET(param_cMOSFET *p, param_solver *ps)
{
  p->radius=get_radius(FET_params);
  p->Lg=FET_params.Lg;
  p->eps_semi=FET_params.eps_s;
  p->Rs=params_parasitic.Rs;
  p->Rd=params_parasitic.Rd;
  p->Cox=get_Cox(FET_params);
  p->temp=temperature;
  p->ni=FET_params.subclass.params_cMOSFET_data->ni;
  p->dphi=FET_params.subclass.params_cMOSFET_data->dphi;
  p->tox=FET_params.tox;
  p->eps_ox=FET_params.eps_ox;
  p->mue=FET_params.mobility;

  ps->left = solver_params.left;
  ps->right = solver_params.right;
}

void get_global_cMESFET(param_cMESFET *p, param_solver *ps)
{
  p->temp=temperature;
  p->radius=get_radius(FET_params);
  p->Lg=FET_params.Lg;
  p->eps_semi=FET_params.eps_s;
  p->Rs=params_parasitic.Rs;
  p->Rd=params_parasitic.Rd;
  p->mue=FET_params.mobility;
  /* p->Cox=FET_params.Cox; */
  /* p->ni=FET_params.ni; */
  /* p->dphi=dphi; */
  /* p->tox=FET_params.tox; */
  /* p->eps_ox=FET_params.eps_ox; */
  p->Nd=FET_params.subclass.params_cMESFET_data->Nd;
  p->Vbi=FET_params.subclass.params_cMESFET_data->Vbi;
  
  ps->left = solver_params.left;
  ps->right = solver_params.right;
}

// cMOSFET: Charges

number frf_Q_cMOSFET(number qq,number V, number Vgs,params_NWFET FET_params)
{
  param_cMOSFET p;
  p.radius=get_radius(FET_params);
  p.Lg=FET_params.Lg;
  p.eps_semi=FET_params.eps_s;
  p.Rs=params_parasitic.Rs;
  p.Rd=params_parasitic.Rd;
  p.Cox=get_Cox(FET_params);
  p.temp=temperature;
  p.ni=FET_params.subclass.params_cMOSFET_data->ni;
  p.dphi=FET_params.subclass.params_cMOSFET_data->dphi;
  p.tox=FET_params.tox;
  p.eps_ox=FET_params.eps_ox;
  p.mue=FET_params.mobility;
  return(func_rootfind_Q_cMOSFET(qq,V, Vgs, p));
}
						  
number Qcharge_cMOSFET(number Vgs)
{
  param_cMOSFET p;
  param_solver ps;
  get_global_cMOSFET(&p, &ps);
  return(func_Qcharge_cMOSFET(Vgs,0,p,ps));
xb}

number Qcharge2_cMOSFET(number Vgs)
{
  param_cMOSFET p;
  param_solver ps;
  get_global_cMOSFET(&p,&ps);
  return(func_Qcharge2_cMOSFET(Vgs,0,p));
}

// cMOSFET: current

#define EXPORT_CMOSFET(name)						 \
  number name(number Vds, number Vgs)				 \
  {													 \
	param_cMOSFET p;								 \
	param_solver ps;								 \
	get_global_cMOSFET(&p,&ps);						 \
	return( func_##name(Vds,Vgs,p,ps));				 \
  }

EXPORT_CMOSFET(Ids_cMOSFET);
EXPORT_CMOSFET(Ids_cMOSFET_R);

#define EXPORT2_CMOSFET(name)						 \
  number name(number Vds, number Vgs)				 \
  {													 \
	param_cMOSFET p;								 \
	param_solver ps;								 \
	get_global_cMOSFET(&p,&ps);						 \
	return( func_##name(Vds,Vgs,p));				 \
  }

EXPORT2_CMOSFET(Ids2_cMOSFET);
EXPORT2_CMOSFET(Ids2_cMOSFET_R);

// cMESFET: current

#define EXPORT_CMESFET(name)				 \
  number name(number Vds, number Vgs)		 \
  {											 \
	param_cMESFET p;						 \
	param_solver ps;						 \
	get_global_cMESFET(&p,&ps);				 \
	return(func_##name(Vds,Vgs,p,ps));		 \
  }

EXPORT_CMESFET(Ids_cMESFET);
EXPORT_CMESFET(Ids_cMESFET_R);

/* number Ids_cMOSFET(number Vds,number Vgs) */
/* { */
/*   param_cMOSFET p; */
/*   get_global_cMOSFET(&p); */
/*   return(Ids0_cMOSFET(Vds,Vgs,p)); */
/* } */

/* number Ids_cMOSFET_R(number Vds,number Vgs) */
/* { */
/*   param_cMOSFET p; */
/*   get_global_cMOSFET(&p); */
/*   return(Ids0_cMOSFET_R(Vds,Vgs,p)); */
/* } */

/* number Ids_cMESFET(number Vds,number Vgs) */
/* { */
/*   param_cMESFET p; */
/*   get_global_cMESFET(&p); */
/*   return(Ids0_cMESFET(Vds,Vgs,p)); */
/* } */

/* number Ids_cMESFET_R(number Vds,number Vgs) */
/* { */
/*   param_cMESFET p; */
/*   get_global_cMESFET(&p); */
/*   return(Ids0_cMESFET_R(Vds,Vgs,p)); */
/* } */
