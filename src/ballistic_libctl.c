/*
 *  ballistic_libctl.c - Time-stamp: <Tue Jul 30 14:55:58 JST 2019>
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
 *  $Id: ballistic_libctl.c 2019-07-30 14:23:51 jmotohisa $
 */

/*! 
  @file ballistic_libctl.c 
  @brief 
  @author J. Motohisa
  @date
*/

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <tgmath.h>
#include "ctl-io.h"

#include "density1d.h"
#include "ballistic.h"

#define GLOBAL_VALUE_DEFINE
// #include "ballistic_libctl.h"

/*!
  @brief
  @param[in]
  @param[out]
  @param[in,out]
  @return
*/

number E0_func(params_ballisticFET0_type pp)
{
  param_ballistic p;
  p.EFermi = pp.Fermi_Energy;
  p.VDS = pp.Vds;
  p.VGS = pp.Vgs;
  p.alpha_D = pp.alpha_D;
  p.alpha_G = pp.alpha_G;
  p.Ceff = pp.C_eff;
  p.alpha = pp.alpha_NP;
  p.ems = pp.effective_mass;
  p.W1 = pp.size_W1;
  p.temp = temperature;
  p.W2 = pp.size_W2;
  p.nmax = pp.n_max;
  p.mmax = pp.m_max;
  return(find_E0(p));
}
  
number funcval_E00(double ene0,params_ballisticFET0_type pp)
{
  param_E0 p;
  param_density1d_all p_density1d_all;

  p_density1d_all.alpha=pp.alpha_NP;
  p_density1d_all.ems=pp.effective_mass;
  p_density1d_all.temp=temperature;
  p_density1d_all.W1=pp.size_W1;
  p_density1d_all.W2=pp.size_W2;
  p_density1d_all.nmax=pp.n_max;
  p_density1d_all.mmax=pp.m_max;

  p.p_density1d_all = p_density1d_all;
  p.EFermi = pp.Fermi_Energy;
  p.VDS = pp.Vds;
  p.VGS = pp.Vgs;
  p.alpha_D = pp.alpha_D;
  p.alpha_G = pp.alpha_G;
  p.Ceff = pp.C_eff;
  return(ene0-func_E00(ene0,&p));
}
