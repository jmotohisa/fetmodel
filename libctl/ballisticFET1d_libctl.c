/*
 *  ballisticFET1d_libctl.c - Time-stamp: <Wed Aug 21 08:53:17 JST 2019>
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
 *  $Id: ballisticFET1d_libctl.c 2019-08-19 20:12:52 jmotohisa $
 */

/*! 
  @file ballisticFET1d_libctl.c 
  @brief 
  @author J. Motohisa
  @date
*/

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <tgmath.h>
#include <gsl/gsl_const_mks.h>
#include <ctl.h>
#include "ctl-io.h"

#include "../src/density1d.h"
#include "../src/ballistic.h"
#include "../src/ballistic_common.h"
#include "cfet_libctl.h"

/* #define GLOBAL_VALUE_DEFINE */
/* #include "ballisticFET1d_libctl.h" */

/*!
  @brief
  @param[in]
  @param[out]
  @param[in,out]
  @return
*/


/*!
  @brief calculate electron density in 1D quantum wire
  @param[in] EFermi : Fermi energy
  @param[in] NW_params
  @return electron density
*/

number density1d(number EFermi, params_NWFET NW_params)
{
  number ems = NW_params.effective_mass;
  //  number Eg = NW_params.bandgap;
  number alphaNP = NW_params.alpha_nonparabolicity;
  integer nmax=NW_params.n_max;
  integer mmax=NW_params.m_max;
  number W1 = get_radius(NW_params);
  number W2 = get_radius2(NW_params);
  double Enm;
  number n1d;
  int n,m;

  switch(NW_params.sizes.which_subclass) {
  case NW_RADIAL:
	if(NW_params.nonparabolicp)
	  {
		n1d=0;
		for(n=1;n<=nmax;n++) {
		  Enm = Ep_n_radial1d(ems,W1,n);
		  n1d += density1d0(EFermi,Enm,ems,temperature);
		}
	  }
	else
	  {
		n1d = 0;
		double Enm0,gamma_nm,alpha_nm,ems_nm;
		for(n=1;n<=nmax;n++) {
		  Enm0 = Ep_n_radial1d(ems,W1,n);
		  gamma_nm = gamma_nm_NP(Enm0,alphaNP);
		  Enm = E_nm_NP(alphaNP,gamma_nm);
		  alpha_nm = alpha_nm_NP(alphaNP,gamma_nm);
		  ems_nm = ems_nm_NP(ems,gamma_nm);
		  n1d += density1d_NP0(EFermi,Enm,alpha_nm,ems_nm,temperature);
		}
	  }
	break;
  case NW_RECT:
  default:
	if(NW_params.nonparabolicp)
	  {
		n1d=0;
		for(n=1;n<=nmax;n++)
		  for(m=1;m<mmax;m++) {
			Enm = Ep_nm_rect1d(ems,W1,W2,n,m);
			n1d += density1d0(EFermi,Enm,ems,temperature);
		  }
	  }
	else
	  {
		n1d = density1d_rect1d_all0(EFermi,ems,temperature,W1,W2,nmax,mmax);
	  }
	break;
  }
  return(n1d);
}

number E0_rect1d(number VDS, number VGS)
{
  double E0;
  params_NWFET p=FET_params ;
  params_ballisticFET b=ballistic_params;
  double W1 = get_radius(p);
  double W2 = get_radius2(p);
  double Cox = Cox_rect(p.eps_ox,p.tox,W1,W2);
  double Cc  = Cc_rect(p.eps_s,W1,W2);
  double Ceff = Cox*Cc/(Cox+Cc);
  // b.C_eff is overridden
  
  E0=E0_rect1d_root0(b.Fermi_Energy, VDS, VGS,
					 b.alpha_D, b.alpha_G, Ceff,
					 p.alpha_nonparabolicity, p.effective_mass, temperature,
					 W1, W2, p.n_max, p.m_max);
  return(E0);
}				 

number Ids_ballistic_rect1d(number VDS, number VGS, number EFs)
{
  params_NWFET p=FET_params ;
  params_ballisticFET b=ballistic_params;
  double W1 = get_radius(p);
  double W2 = get_radius2(p);
  double Cox = Cox_rect(p.eps_ox,p.tox,W1,W2);
  double Cc  = Cc_rect(p.eps_s,W1,W2);
  double Ceff = Cox*Cc/(Cox+Cc);
  // b.C_eff is overridden
  double ids;
  ids=Ids_ballistic1d_rect1dNP0(VDS, VGS, EFs,
								b.Fermi_Energy,
								b.alpha_D, b.alpha_G, Ceff,
								p.alpha_nonparabolicity, p.effective_mass, temperature,
								W1, W2, p.n_max, p.m_max);
  return(ids);
}
