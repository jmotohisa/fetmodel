/*
 *  qw2d_test.c - Time-stamp: <Thu Aug 31 09:25:54 JST 2023>
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
 *  $Id: qw2d_test.c 2023-08-31 07:47:01 jmotohisa $
 */

/*! 
  @file qw2d_test.c 
  @brief 
  @author J. Motohisa
  @date
*/

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <tgmath.h>

/* #define GLOBAL_VALUE_DEFINE */
/* #include "qw2d_test.h" */

/*!
  @brief
  @param[in]
  @param[out]
  @param[in,out]
  @return
*/

#include "../src/density1d.h"
#include "../src/density2d.h"
#include "../src/ballistic.h"

  /* GLOBAL double Ids_ballistic2d_QW0(double VDS, double VGS, double EFs,  */
  /* 				    double alpha_D, double alpha_G, */
  /* 				    double Ceff, */
  /* 				    double alpha, double ems, double temp, */
  /* 				    double W1, int nmax); */

int main()
{
  double VDS, VGS, EFs;
  double alpha_D, alpha_G;
  double Ceff;
  double alpha, ems,  temp;
  double W1;
  int nmax;

  double Eg = 0.36;
  double epsOX = 8.5;
  double epsS = 8.9;
  double tOX = 20e-9;
  
  temp = 300;
  ems = 0.067;
  W1 = 5*1e-9;
  double W2 = 8e-9;//  # does not matter but set value for clarity
  alpha = alpha_NP(Eg, ems);
  double Cox = epsOX*8.85e-12/tOX;
  double Cc = sqrt(1.6e-19*1e21/(2*epsS*8.86e-12*1));
  alpha_D = 0;
  alpha_G = 1;
  nmax=10;
  EFs=0;
  
  /* Vgsmin=-0.1; */
  /* Vgsmax=0.1; */
  /* dVgs=0.01 */
  double ids;
  VDS=0.5;
  VGS=-0.1;
  ids=Ids_ballistic2d_QW0(VDS, VGS, EFs,
			  alpha_D, alpha_G,
			  Cox,
			  alpha, ems, temp,
			  W1, nmax);
  printf("%le\n",ids);
  return 1;
}
