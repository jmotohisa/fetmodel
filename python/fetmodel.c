/*
 *  fetmodel.c - Time-stamp: <Tue Mar 10 17:45:14 JST 2020>
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
 *  $Id: fetmodel.c 2019-07-12 14:25:01 jmotohisa $
 */

/*! 
  @file fetmodel.c 
  @brief 
  @author J. Motohisa
  @date
*/

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <tgmath.h>
#include "../src/ccm.h"
#include "../src/density1d.h"
#include "../src/ballistic.h"

#define GLOBAL_VALUE_DEFINE
#include "fetmodel.h"

/*!
  @brief
  @param[in]
  @param[out]
  @param[in,out]
  @return
*/

/* #define DEFINE_CONSTRUCTOR(STRUCT)				\ */
/*   STRUCT * STRUCT_new()						\ */
/*   {												\ */
/* 	STRUCT *p = malloc(sizeof(STRUCT));			\ */
/* 	if(p ==NULL) {								\ */
/* 	  perror("malloc");							\ */
/* 	  exit(1);									\ */
/* 	}											\ */
/* 	return(p);									\ */
/*   } */

/* param_cMOSFET *param_cMOSFET_new() */
/* { */
/*   param_cMOSFET *p = malloc(sizeof(param_cMOSFET)); */
/*   if(p ==NULL) { */
/* 	perror("malloc"); */
/* 	exit(1); */
/*   } */
/*   return(p); */
/* } */

/* param_cMESFET *param_cMESFET_new() */
/* { */
/*   param_cMESFET *p = malloc(sizeof(param_cMESFET)); */
/*   if(p ==NULL) { */
/* 	perror("malloc"); */
/* 	exit(1); */
/*   } */
/*   return(p); */
/* } */

// Cylindrical MOSFET
// charges:: numpy compatible
void Qapprox_cMOS_func(double *in_array,double *out_array,int size,
					   param_cMOSFET p)
{
  int i;
  for(i=0;i<size;i++)
	*(out_array+i)=func_Qcharge2_cMOSFET(0,*(in_array+i),p);
}

void Q_cMOS_func(double *in_array,double *out_array,int size,
				 param_cMOSFET p)
{
  int i;
  param_solver ps;
  for(i=0;i<size;i++)
	*(out_array+i)=func_Qcharge_cMOSFET(0,*(in_array+i),p,ps);
}

// Current:: numpy Vgs
void Ids0_cMOS_func(double *in_array, double *out_array, int size,
				   double Vds,param_cMOSFET p)
{
  int i;
  for(i=0;i<size;i++)
	*(out_array+i)=func_Ids2_cMOSFET(Vds,*(in_array+i),p);
}

void Ids_cMOS_func(double *in_array, double *out_array, int size,
				   double Vds,param_cMOSFET p)
{
  int i;
  param_solver ps;
  for(i=0;i<size;i++)
	*(out_array+i)=func_Ids_cMOSFET(Vds,*(in_array+i),p,ps);
}

void Ids0_cMOS_R_func(double *in_array, double *out_array, int size,
				   double Vds,param_cMOSFET p)
{
  int i;
  for(i=0;i<size;i++)
	*(out_array+i)=func_Ids2_cMOSFET_R(Vds,*(in_array+i),p);
}

void Ids_cMOS_R_func(double *in_array, double *out_array, int size,
				   double Vds,param_cMOSFET p)
{
  int i;
  param_solver ps;
  for(i=0;i<size;i++)
	*(out_array+i)=func_Ids_cMOSFET_R(Vds,*(in_array+i),p,ps);
}

// cyrlindircal MESFET
void Ids_cMES_func(double *in_array, double *out_array, int size,
				   double Vds,param_cMESFET p)
{
  int i;
  param_solver ps;
  for(i=0;i<size;i++)
	*(out_array+i)=func_Ids_cMESFET(Vds,*(in_array+i),p,ps);
}

void Ids_cMES_R_func(double *in_array, double *out_array, int size,
					 double Vds,param_cMESFET p)
{
  int i;
  param_solver ps;
  for(i=0;i<size;i++)
	*(out_array+i)=func_Ids_cMESFET_R(Vds,*(in_array+i),p,ps);
}

// ballistic FET
/* // DEFINE_CONSTRUCTOR(param_ballistic); */
/* param_ballistic *param_ballistic_new() */
/* { */
/*   param_ballistic *p = malloc(sizeof(param_ballistic)); */
/*   if(p ==NULL) { */
/* 	perror("malloc"); */
/* 	exit(1); */
/*   } */
/*   return(p); */
/* } */
