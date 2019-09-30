/*
 *  capacitor.c - Time-stamp: <Mon Aug 19 12:31:17 JST 2019>
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
 *  $Id: capacitor.c 2019-08-10 10:10:14 jmotohisa $
 */

/*! 
  @file capacitor.c 
  @brief Capacitors in various structures
  @author J. Motohisa
  @date
*/

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <tgmath.h>
#include <gsl/gsl_const_mks.h>

#include "ballistic_common.h"

#define GLOBAL_VALUE_DEFINE
#include "capacitor.h"

/*!
  @brief
  @param[in]
  @param[out]
  @param[in,out]
  @return
*/

double Cox_rect(double epsOX,double tOX,double W1, double W2)
{
  return(2*EPSILON(epsOX)*(W1+W2)/tOX + 2.232*EPSILON(epsOX));
}

double Cox_rect_area(double epsOX,double tOX,double W1, double W2)
{
  return(Cox_rect(epsOX,tOX,W1,W2)/(2*(W1+W2)));
}

double Cc_rect(double epsS, double W1, double W2)
{
  return (6.94*EPSILON(epsS)*(W1*W1+W2*W2)/(W1*W2));
}

double Cox_radial(double epsOX, double tOX, double radius)
{
  return(2*M_PI*EPSILON(epsOX)/(log(1+tOX/radius)));
}

double Cox_radial_area(double epsOX, double tOX, double radius)
{
  return(EPSILON(epsOX)/(radius*log(1+tOX/radius)));
}

