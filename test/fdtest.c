/*
 *  fdtest.c - Time-stamp: <Thu Aug 31 09:07:05 JST 2023>
 *
 *   Copyright (c) 2023  jmotohisa (Junichi Motohisa)  <motohisa@ist.hokudai.ac.jp>
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
 *  $Id: fdtest.c 2023-08-31 08:47:09 jmotohisa $
 */

/*! 
  @file fdtest.c 
  @brief 
  @author J. Motohisa
  @date
*/

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <tgmath.h>

/* #define GLOBAL_VALUE_DEFINE */
/* #include "fdtest.h" */

/*!
  @brief
  @param[in]
  @param[out]
  @param[in,out]
  @return
*/

#include <gsl/gsl_sf_fermi_dirac.h>
#include "../src/fd.h"

int main()
{
  double ene;
  double fd01,fd02,fd03;
  double fdh1,fdh2,fdh3;
  double fdm1,fdm2,fdm3;
  
  ene=-100;
  fd01=gsl_sf_fermi_dirac_0(ene);
  fd02=fd_0(ene);
  fd03=fd01;
  
  fdh1=gsl_sf_fermi_dirac_half(ene);
  fdh2=fd_half(ene);
  fdh3=exp(ene);
  
  fdm1=fd_mhalf(ene);
  fdm2=gsl_sf_fermi_dirac_mhalf(ene);
  fdm3=exp(ene);

  /* printf("%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\n", */
  /* 		 ene,fd01,fd02,fd03,fdh1,fdh2,fdh3,fdm1,fdm2,fdm3); */
  printf("%le\t%le\t%le\t%le\t%le\t%le\n",ene,fdh1,fd01,fdm1,fdh1/fd01,fd01/fdm1);
}
