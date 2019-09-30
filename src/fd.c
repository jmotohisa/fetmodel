/*
 *  fd.c - Time-stamp: <Mon Sep 16 22:17:49 JST 2019>
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
 *  $Id: fd.c 2019-09-16 21:45:41 jmotohisa $
 */

/*! 
  @file fd.c 
  @brief 
  @author J. Motohisa
  @date
*/

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <tgmath.h>
#include <gsl/gsl_sf_fermi_dirac.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_machine.h> 
#define GLOBAL_VALUE_DEFINE
#include "fd.h"

/*!
  @brief
  @param[in]
  @param[out]
  @param[in,out]
  @return
*/

double fd_mhalf(double x)
{
  gsl_sf_result result;
  int status;
  double retval;
  if(x < GSL_LOG_DBL_MIN)
	return 0;
  status=gsl_sf_fermi_dirac_mhalf_e(x,&result);
  /* val=gsl_sf_fermi_dirac_mhalf (x); */
  retval=result.val;
  if(status!=GSL_SUCCESS) {
	printf ("Error in gsl_sf_fermi_dirac_mhalf: %s\n", gsl_strerror (status));
	switch (status) {
	case GSL_EUNDRFLW:
	  retval=0.;
	  break;
	case GSL_EOVRFLW:
	default:
	  retval=2*sqrt(x/M_PI);
	  break;
	}
  }
	return(retval);
}

double fd_0(double x)
{
  gsl_sf_result result;
  int status;
  double retval;
  status=gsl_sf_fermi_dirac_0_e(x,&result);
  /* val=gsl_sf_fermi_dirac_0 (x); */
  retval=result.val;
  if(status!=GSL_SUCCESS) {
	printf ("Error in gsl_sf_fermi_dirac_0: %s\n", gsl_strerror (status));
	switch (status) {
	case GSL_EUNDRFLW:
	  retval=0.;
	  break;
	case GSL_EOVRFLW:
	default:
	  retval=x;
	  break;
	}
  }
	return(retval);
}


double fd_half(double x)
{
  gsl_sf_result result;
  int status;
  double retval;
  status=gsl_sf_fermi_dirac_half_e(x,&result);
  /* val=gsl_sf_fermi_dirac_half (x); */
  retval=result.val;
  if(status!=GSL_SUCCESS) {
	printf ("Error in gsl_sf_fermi_dirac_half: %s\n", gsl_strerror (status));
	switch (status) {
	case GSL_EUNDRFLW:
	  retval=0.;
	  break;
	case GSL_EOVRFLW:
	default:
	  retval=4.*sqrt(x/M_PI)*x/3.;
	  break;
	}
  }
	return(retval);
}
