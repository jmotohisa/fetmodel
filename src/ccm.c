/*
 *  ccm.c - Time-stamp: <Sun Mar 08 19:48:52 JST 2020>
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
 *  $Id: ccm.c 2019-07-12 16:04:01 jmotohisa $
 */

/*! 
  @file ccm.c 
  @brief FET models for Cylindrical MOSFET and Cylindrical MESFET
  @author J. Motohisa
  @date
*/

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <tgmath.h>

#define GLOBAL_VALUE_DEFINE
#include "ccm.h"

/*!
  @brief 
  @param[in]
  @param[out]
  @param[in,out]
  @return
*/

// 
// in Gradual Channel Approximation

/*
  MOSFET: Iniguez, B., Jimenez, D., Roig, J., Hamid, H., Marsal, L., & Pallares, J. (2005).
   Explicit continuous model for long-channel undoped surrounding gate MOSFETs.
   IEEE Transactions on Electron Devices, 52(8), 1868–1873.

  MESFET: see research report JM
 */

#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#include <math.h>
#include "PhysicalConstants.h"

double Ids00_cMOSFET(double , double , param_cMOSFET);
double n_intrinsic(double , double , double ,double ,double );
double gm0_cMESFET(double ,double , param_cMESFET, param_solver );
double r2root0(double ,double , param_cMESFET , param_solver );
double Q_approx0(double , double , double , double , param_cMOSFET );

#define POW2(X) (X)*(X)
#define kBT (KBC*p.temp)
#define Vth (kBT/EC)
#define Q0 (4*p.eps_semi*EPSILON/p.radius*Vth)
#define delta (EC*EC*p.ni/(kBT*p.eps_semi*EPSILON))
#define V0 (p.dphi + Vth*log(8/(delta*p.radius*p.radius)))
#define KcMOSFET (2*PI * p.radius/p.Lg *p.mue)
#define KcMESFET (PI*POW2(EC*p.Nd)*p.mue*POW2(POW2(p.radius))/(16*EPSILON*p.eps_semi*p.Lg))

#define Vth_cMESFET (p.Vbi-EC*p.Nd*POW2(p.radius)/(4*EPSILON*p.eps_semi))
#define  K2cMESFET (PI*EC*p.Nd*p.mue*POW2(p.radius)/(p.Lg))

// Struct types for the solution of equation with GSL
struct func_Qcharge_cMOSFET_param {double V ; double Vgs; param_cMOSFET p;};
struct R_cMESFET_param {double V ; double Vgs; param_cMESFET p;};
struct Ids_params {double Vds; double Vgs;param_cMOSFET cMOS;param_cMESFET cMES;param_solver ps;};

// solution for cylindrical MOSFET
double qroot_newton(double, double, param_cMOSFET, param_solver );
double qroot_brent(double, double, param_cMOSFET, param_solver );
double qfunc_cMOSFET_gsl(double, void *);
double qdfunc_cMOSFET_gsl(double, void *);
void qfdfunc_cMOSFET_gsl(double, void *, double *, double *);
double r_cMESFET_gsl(double, void *);

int Ids0_cMOSFET_RmodFunc(gsl_vector *, void *, gsl_vector *);
int Ids_cMOSFET_RmodFunc(gsl_vector *, void *, gsl_vector *);
int Ids_cMESFET_RmodFunc(gsl_vector *, void *, gsl_vector *);

/* void recalc_param_cMOSFET() */
/* { */
/*   // temp, eps_semi, dphi,Cox and ni should be calculated before calling this proc */
/*   //  Vth=temp*KBC/EC; */
/*   //	Cox = eps_ox*EPSILON/(radius*log(1+tox/radius)) */
/*   //  Q0=4*eps_semi*EPSILON/radius*Vth; */
/*   //	ni = n_intrinsic(temperature,me,mh,Eg,6) */
/*   //  delta=(EC*ni/(KBT*(eps_semi*EPSILON))); */
/*   //  V0= (dphi + Vth*log(8/(delta*POW2(radius)))); */
/*   //  K = 2*PI * radius/Lg *mue; */
/* } */

/* void SetParam_cMOSFET() */
/* {	  */
/*   //  double radius=50,Lg=1,Cox=400,mue=400,ni=1.45e10,temp=300,eps_semi=11.9,dphi=0; */
/*   /\* */
/* 	Prompt radius,"nanowire radius (nm)" */
/* 	Prompt Lg,"Gate Legth (micron)" */
/* 	Prompt Cox, "Oxide capacitance (pF/m)" */
/* 	Prompt mue,"mobility (cm^2/Vs)" */
/* 	Prompt ni,"intrinsic carrier concentration (cm^-3)" */
/* 	Prompt temp,"temperature (K)" */
/* 	Prompt eps_semi,"dielectric constant" */
/* 	Prompt dphi,"work function difference" */
/* 	PauseUpdate;Silent 1 */
/*   *\/ */

/*   //  temp=temp; */
/*   radius=radius*1e-9; */
/*   Lg=Lg*1e-6; */
/*   Cox=Cox*1e-12; */
/*   mue=mue*1e-4; */
/*   ni=ni*1e6; */
/*   //  eps_semi=eps_semi; */
/*   //  dphi = dphi; */
/*   recalc_param_cMOSFET(); */
/* } */

/* // intrinsic carrier concentration */
/* void SetParam1_cMESFET(int flg) */
/* { */
/*   //  double flg=1,ni=1.48e10,temp=300,me,mh,Eg,gv=6; */
/*   /\* */
/* 	Prompt flg,"Calculate ni",popup,"yes;no" */
/* 	Prompt ni,"intrinsic carrier concentration (cm^-3)" */
/* 	Prompt temp,"temperature (K)" */
/* 	Prompt me,"DOS effective mass of CB (m0)" */
/* 	Prompt me,"DOS effective mass of VB (m0)" */
/* 	Prompt Eg,"Bandgap (eV)" */
/* 	Prompt gv,"valley degeneracy (CB)" */
/*   *\/ */
	
/* 	if(flg==1) */
/* 	  ni=n_intrinsic(temp,me_dos,mh_dos,Eg,gv_e); */
/* 	else */
/* 	  ni=ni*1e6; */

/* 	ni=ni; */
/* 	recalc_param_cMOSFET(); */
/* } */

// gate capacitance
void SetParam2_cMOSFET(int flg, param_cMOSFET *p)
{
  //  double flg=1,Cox=400,radius=50,tox=1.5,eps_ox=3.9;
	  /*	Prompt flg,"Calculate Cox",popup,"yes;no"
	Prompt Cox,"gate capacitance (pF/m)"
	Prompt radius,"nanowire radius (nm)"
	Prompt tox,"oxide thickness (nm)"
	Prompt eps_ox,"dielectric const of oxide"
	  */
	if(flg==1)
	  p->Cox = p->eps_ox*EPSILON/(p->radius*log(1+p->tox/p->radius));
	else
	  p->Cox = p->Cox*1e-12;

	/* recalc_param_cMOSFET(); */
}

/* void param_Si_cMOSFET() */
/* { */
/*   //  me_dos = (0.98*POW2(0.19))^(1/3); */
/*   //  mh_dos=(0.168^(3/2)+0.498^(3/2))^(2/3); */
/*   Eg=1.12; */
/*   eps_semi = 11.9; */
/*   ni = 1.45e10 *1e6; */
/*   gv_e=6; */
/* //	Nc = 2.8e19 *1e6; */
/* //	Nv = 1.04e19 *1e6; */
/* } */

/* void param_Ge_cMOSFET() */
/* { */
/*   //  me_dos = (0.98*POW2(0.19))^(1/3); */
/*   //  mh_dos=(0.168^(3/2)+0.498^(3/2))^(2/3); */
/*   Eg=1.12; */
/*   eps_semi = 11.9; */
/*   ni = 1.45e10 *1e6; */
/*   gv_e=8; */
/* //	me = (1.64*POW2(0.082))^(1/3) */
/* //	mh=(0.48^(3/2)+0.288^(3/2))^(2/3) */
/* //	Eg=0.66 */
/* //	eps_semi = 16.0 */
/* //	ni = 2.4e13 *1e6 */
/* //	Nc = 1.04e19 *1e6 */
/* //	Nv = 6.0e18 *1e6 */
/* } */

/* void param_GaAs_cMOSFET() */
/* { */
/*   //  me_dos = (0.98*POW2(0.19))^(1/3); */
/*   //  mh_dos=(0.168^(3/2)+0.498^(3/2))^(2/3); */
/*   Eg=1.12; */
/*   eps_semi = 11.9; */
/*   ni = 1.45e10 *1e6; */
/*   gv_e=1; */
/* //	me = 0.067 */
/* //	mh=(0.0828^(3/2)+0.458^(3/2))^(2/3) */
/* //	Eg=1.43 */
/* //	eps_semi = 13.1 */
/* //	ni = 1.79e16 *1e6 //@300K */
/* //	Nc = 4.7e17 *1e6 */
/* //	Nv = 7.0e18 *1e6 */
/* } */

/***************************************************/
/***************************************************/

/*!
  @brief Approximate form of Q (eq.17 in Iniguez Trans. ED)
  @param[in] V: Drain voltage
  @param[in] Vgs : Gate voltage
  @param[in] p: Parameters for cMOSFET
  @return
*/
double Q_approx(double V, double Vgs, param_cMOSFET p)
{
  double qp0,Vt0,dVt0;
  qp0 = Q_approx0(V,Vgs,V0,0,p);
  Vt0=V0 + 2*Vth*log(1+qp0/Q0);
  dVt0 = (2*p.Cox*POW2(Vth)/Q0)*qp0/(Q0+qp0);
  return(Q_approx0(V,Vgs,Vt0,dVt0,p));
}

/*!
  @brief Unified approximated formula for Q (eq. 13 in Iniguez et al., Trans. ED)
  @param[in] V, Vgs, Vt, deltaVt
  @param[in] p: Parameters for cMOSFET
  @return
*/
double Q_approx0(double V,double Vgs,double Vt,double deltaVt, param_cMOSFET p)
{
  double a,b,c1,c2;
  a=2*p.Cox*Vth*Vth/Q0;
  b=2*Vth*log(1+exp((Vgs-Vt+deltaVt-V)/(2*Vth)));
  c1=sqrt(a*a+b*b)-a;
  c2=b*b/(2*a);
  if(c1 < DBL_EPSILON)
	return(p.Cox*c2);
  else
	return(p.Cox*c1);
}

//find solution for Q usint GSL library
double qroot0(double V, double Vgs, param_cMOSFET p, param_solver ps)
{
  /* return(qroot_brent1(V,Vgs,p,ps)); */
  /* return(qroot_brent2(V,Vgs,p)); */
  return(qroot_newton(V,Vgs,p,ps));
}

// Solver with Blent:
double qroot_brent0(double V,double Vgs,param_cMOSFET p,double low, double high, int max_iter)
{
  gsl_function F;
  struct func_Qcharge_cMOSFET_param params={V,Vgs,p};
  int status;
  int iter = 0;
  double r;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  
  F.function = &qfunc_cMOSFET_gsl;
  F.params = &params;
  //FindRoots/Q/L=(low) qfunc_cMOSFET,param_cMOSFET;

  T= gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc(T);

  gsl_root_fsolver_set(s,&F,low,high);
  do {
	iter++;
	status = gsl_root_fsolver_iterate(s);
	r=gsl_root_fsolver_root(s);
	low = gsl_root_fsolver_x_lower(s);
	high = gsl_root_fsolver_x_upper(s);
	status = gsl_root_test_interval(low,high,0,0.001);
  } while (status==GSL_CONTINUE && iter < max_iter);
  
  return(r);
}
double qroot_brent1(double V,double Vgs,param_cMOSFET p,param_solver ps)
{
  int max_iter = 100;
  double low=ps.left,high=ps.right;
  return(qroot_brent0(V,Vgs,p,low,high,max_iter));
}

double qroot_brent2(double V,double Vgs,param_cMOSFET p)
{
  int max_iter = 100;
  double low,high;
  double qqq;
  qqq=Q_approx(V,Vgs,p);
  low=qqq/100.;
  high=qqq*100.;
  return(qroot_brent0(V,Vgs,p,low,high,max_iter));
}

  
// Solver with newton:
double qroot_newton(double V,double Vgs, param_cMOSFET p, param_solver ps)
{
  int status;
  int iter=0,max_iter=100;
  const gsl_root_fdfsolver_type *T;
  gsl_root_fdfsolver *s;
  double x0,x;
  gsl_function_fdf FDF;
  struct func_Qcharge_cMOSFET_param params; //={V,Vgs,p};
  params.V=V;
  params.Vgs=Vgs;
  params.p=p;
  
  x=Q_approx(V,Vgs,p);

  FDF.f = &qfunc_cMOSFET_gsl;
  FDF.df = &qdfunc_cMOSFET_gsl;
  FDF.fdf = &qfdfunc_cMOSFET_gsl;
  FDF.params = &params;
  T = gsl_root_fdfsolver_newton;
  s = gsl_root_fdfsolver_alloc(T);
  gsl_root_fdfsolver_set(s,&FDF,x);

  do {
	iter++;
	status = gsl_root_fdfsolver_iterate(s);
	x0=x;
	x=gsl_root_fdfsolver_root(s);
	status = gsl_root_test_delta(x,x0,0,1e-3);
  } while (status==GSL_CONTINUE && iter < max_iter);
  
  return(x);
}

double qfunc_cMOSFET(double qq,double V, double Vgs, param_cMOSFET p)
{
  double qqq1,qqq2;
  // dirty hack
  if(qq<0)
	qq=fabs(qq)/10;
  qqq1 = Vgs-p.dphi-V-Vth*log(8/(delta*POW2(p.radius)));
  qqq2 =(qq/p.Cox + Vth*(log(qq/Q0)+log(1+(qq/Q0))));
  return(qqq1-qqq2);
}

double qfunc_cMOSFET_gsl(double qq, void *pp)
{
  struct func_Qcharge_cMOSFET_param *params = (struct func_Qcharge_cMOSFET_param *) pp;
  double V = (params->V);
  double Vgs = (params->Vgs);
  param_cMOSFET p = (params->p);
  return(qfunc_cMOSFET(qq,V,Vgs,p));
}

double qdfunc_cMOSFET_gsl(double qq, void *pp)
{
  struct func_Qcharge_cMOSFET_param * params = (struct func_Qcharge_cMOSFET_param *) pp;
  /* double V = (params->V); */
  /* double Vgs = (params->Vgs); */
  param_cMOSFET p = (params->p);
  double qqq1;
  qqq1 =-(1/p.Cox + Vth*(1/qq+1/(qq+Q0)));
  return(qqq1);
}

void qfdfunc_cMOSFET_gsl(double qq, void *pp, double *f, double *df)
{
  struct func_Qcharge_cMOSFET_param * params = (struct func_Qcharge_cMOSFET_param *) pp;
  double V = (params->V);
  double Vgs = (params->Vgs);
  param_cMOSFET p = (params->p);
  *f=qfunc_cMOSFET(qq,V,Vgs,p);
  *df=-(1/p.Cox + Vth*(1/qq+1/(qq+Q0)));
}


// Drain current
/*
  @brief Drain current calculated using Eq.17 in Iniguez Trans. ED
  @param[in] Vds: Drain voltage
  @param[in] Vgs : Gate voltage
  @param[in] p: Parameters for cMOSFET
  @return
*/
double Ids0_cMOSFET(double Vds,double Vgs,param_cMOSFET p)
{
  double QS,QD;
  QS = Q_approx(0,Vgs,p);
  QD = Q_approx(Vds,Vgs,p);
  return(Ids00_cMOSFET(QS,QD,p)*KcMOSFET);
}

double Ids_cMOSFET0(double Vds,double Vgs,param_cMOSFET p,param_solver ps)
{
  double QS,QD;
  QS = qroot0(0,Vgs,p,ps);
  QD = qroot0(Vds,Vgs,p,ps);
  return(Ids00_cMOSFET(QS,QD,p)*KcMOSFET);
}

double Ids00_cMOSFET(double QS,double QD, param_cMOSFET p)
{
  double ids1,ids2;
  ids1 = 2*Vth*(QS-QD)+(POW2(QS)-POW2(QD))/p.Cox;
  ids2 = Vth*Q0*log((QD+Q0)/(QS+Q0));
  return(ids1+ids2);
}

// Effect of source and drain resisitance (1) Approximate solution for Q
double Ids0_cMOSFET_R(double Vds,double Vgs,param_cMOSFET p)
{
  double idsmod;
  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;
  int status;
  /* size_t i; */
  size_t iter=0;
  const size_t n=2;
  gsl_multiroot_function F;
  //  param_cMESFET cMES;
  struct Ids_params params;// = {Vds, Vgs, c,cMES};
  double x_init[2]={1,0};
  gsl_vector *x = gsl_vector_alloc(n);
  params.Vds=Vds;
  params.Vgs=Vgs;
  //  params.ps=ps;
  /* params.cMES=cMES; */
  
  F.f = &Ids0_cMOSFET_RmodFunc;
  F.n = 2;
  F.params = &params;

  gsl_vector_set(x,0,x_init[0]);
  gsl_vector_set(x,1,x_init[1]);
  T=gsl_multiroot_fsolver_hybrids;
  s=gsl_multiroot_fsolver_alloc(T,2);

  gsl_multiroot_fsolver_set(s,&F,x);

  do {
	iter++;
	status = gsl_multiroot_fsolver_iterate(s);
	if(status)
	  break;
	status = gsl_multiroot_test_residual(s->f,1e-7);
  } while(status==GSL_CONTINUE && iter <1000);
  idsmod = Ids0_cMOSFET(gsl_vector_get(s->x,0),gsl_vector_get(s->x,1),p);
  gsl_multiroot_fsolver_free(s);
  gsl_vector_free(x);
  return(idsmod);
}
   
	
int Ids0_cMOSFET_RmodFunc(gsl_vector *x, void *pp, gsl_vector *f)
{
  struct Ids_params * params = (struct Ids_params *)pp;
  const double Vds = (params->Vds);
  const double Vgs = (params->Vgs);
  const param_cMOSFET p = (params->cMOS);
  const double x_Vds = gsl_vector_get(x,0);
  const double x_Vgs = gsl_vector_get(x,1);
  
  gsl_vector_set(f,0,x_Vgs-Vgs+Ids0_cMOSFET(x_Vds,x_Vds,p)*p.Rs);
  gsl_vector_set(f,1,x_Vds-Vds+Ids0_cMOSFET(x_Vds,x_Vds,p)*(p.Rs+p.Rd));

  return GSL_SUCCESS;
}


// effect of source and drain resisitance (2) Rigorous solution for Q
double Ids_cMOSFET_R0(double Vds,double Vgs,param_cMOSFET p,param_solver ps)
{
  double idsmod;
  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;
  int status;
  /* size_t i; */
  size_t iter=0;
  const size_t n=2;
  gsl_multiroot_function F;
  /* param_cMESFET cMES; */
  struct Ids_params params;// = {Vds, Vgs, cMOS,cMES};
  double x_init[2]={1,0};
  gsl_vector *x = gsl_vector_alloc(n);
  params.Vds=Vds;
  params.Vgs=Vgs;
  params.cMOS=p;
  params.ps = ps;
  /* params.cMES=cMES; */
  
  F.f = &Ids_cMOSFET_RmodFunc;
  F.n = 2;
  F.params = &params;

  gsl_vector_set(x,0,x_init[0]);
  gsl_vector_set(x,1,x_init[1]);
  T=gsl_multiroot_fsolver_hybrids;
  s=gsl_multiroot_fsolver_alloc(T,2);

  gsl_multiroot_fsolver_set(s,&F,x);

  do {
	iter++;
	status = gsl_multiroot_fsolver_iterate(s);
	if(status)
	  break;
	status = gsl_multiroot_test_residual(s->f,1e-7);
  } while(status==GSL_CONTINUE && iter <1000);
  idsmod = Ids_cMOSFET(gsl_vector_get(s->x,0),gsl_vector_get(s->x,1),p,ps);
  gsl_multiroot_fsolver_free(s);
  gsl_vector_free(x);
  return(idsmod);
}
   
	
int Ids_cMOSFET_RmodFunc(gsl_vector *x, void *pp, gsl_vector *f)
{
  struct Ids_params * params = (struct Ids_params *)pp;
  const double Vds = (params->Vds);
  const double Vgs = (params->Vgs);
  const param_cMOSFET p = (params->cMOS);
  const param_solver ps = (params->ps);
  const double x_Vds = gsl_vector_get(x,0);
  const double x_Vgs = gsl_vector_get(x,1);
  
  gsl_vector_set(f,0,x_Vgs-Vgs+Ids_cMOSFET(x_Vds,x_Vds,p,ps)*p.Rs);
  gsl_vector_set(f,1,x_Vds-Vds+Ids_cMOSFET(x_Vds,x_Vds,p,ps)*(p.Rs+p.Rd));

  return GSL_SUCCESS;
}


// intrinsic carrier concentration
// n_i = 2 (k_B T \over 2 \PI \hbar^2)^{3/2} (m_e m_h)^{3/4} \exp(-E_g/2 k_B T)
// temp: temperature (K)
// me,mh: denstiy of states effective mass for electron and hole) m_0
// Eg : bandgap (eV)

double n_intrinsic(double temp,double me,double mh,double Eg,double gv)
{
  return(exp(-Eg*EC/(2*KBC*temp)));
  //  return(gv*2*sqrt(KBC*temp*MEL/(2*PI*HBAR*HBAR))^3*(me*mh)^(3/4)*exp(-Eg*EC/(2*KBC*temp)));
}


/////////////////////////
// Cylyndrical MESFET

void SetParams_cMESFET(double radius,double Lg,double Nd,double Vbi,double eps_semi,double mue)
{
  //double radius=50,Lg=1,ND=2.3e18,Vbi=0.34,eps_semi=14.5,mue=400;
	  /*	Prompt radius,"nanowire radius (nm)"
	Prompt Lg,"Gate Legth (micron)"
	Prompt Nd,"DoPIng Density (cm^-3)"
	Prompt Vbi,"Vbi (V)"
	Prompt eps_semi,"Dielectric constant"
	Prompt mue,"mobility (cm^2/Vs)"
	  */

  //  Vbi=Vbi;
  radius=radius* 1e-9;
  //  eps_semi=eps_semi;
  Nd=Nd*1e6;
  mue=mue*1e-4;
  Lg=Lg*1e-6;
  //  KcMESFET=PI*POW2(EC*Nd)*mue*POW2(POW2(radius))/(16*EPSILON*eps_semi*Lg);
  //  Vth=Vbi-EC*Nd*POW2(radius)/(4*EPSILON*eps_semi);
  //  K2=PI*EC*Nd*mue*POW2(radius)/(Lg);
}

double func_Ids_cMESFET(double Vds,double Vgs, param_cMESFET p, param_solver ps)
{
  double r2d,r2s,d0,s0;

  r2d = r2root0(Vds,Vgs,p,ps);
  if(r2d ==0)
	d0=0;
  else
	d0 = POW2(r2d)*(1-2*log(r2d));
  
  if(Vgs>Vth_cMESFET && Vgs<p.Vbi)
	{
	  r2s = r2root0(0,Vgs,p,ps);
	  s0=POW2(r2s)*(1-2*log(r2s));
	}
  else
	{
	  if(Vgs < Vth_cMESFET)
		{		
		  r2s =0;
		  s0=0;
		}
	  else
		{
		  r2s=1;
		  s0=1;
		}
	}
  
  return((s0-d0)*KcMESFET);
}

double gm_cMESFET(double Vds,double Vgs,param_cMESFET p,param_solver ps)
{
  return(K2cMESFET*gm0_cMESFET(Vds,Vgs,p,ps));
}

double gm0_cMESFET(double Vds,double Vgs, param_cMESFET p, param_solver ps)
{
  double r2d,r2s;

  r2d = r2root0(Vds,Vgs, p, ps);
  if(Vgs>Vth_cMESFET && Vgs<p.Vbi)
	r2s = r2root0(0,Vgs,p,ps);
  else
	{
	  if(Vgs < Vth_cMESFET)
		r2s =0;
	  else
		r2s=1;
	  
	}
  return(r2s-r2d);
}

double r2root0(double V,double Vgs, param_cMESFET p, param_solver ps)
{
  double low,high;
  double VbiV,Vbi_r;
  //  double V_root;
  gsl_function F;
  struct R_cMESFET_param params={V,Vgs,p};
  int status;
  int iter = 0, max_iter = 100;
  double r;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;

  F.function = &r_cMESFET_gsl;
  F.params = &params;

  low=ps.left;
  high=ps.right;
  VbiV=p.Vbi-Vgs+V;
  Vbi_r = 4*(p.eps_semi*EPSILON)/(EC*p.Nd*POW2(p.radius))*VbiV;
  if(Vbi_r>0 && Vbi_r<1) // param_cMESFET<1 &&
	{
	  //	  FindRoots/Q/L=(low)/H=(high) r_cMESFET,param_cMESFET;
	  T= gsl_root_fsolver_brent;
	  s = gsl_root_fsolver_alloc(T);
	  
	  gsl_root_fsolver_set(s,&F,low,high);
	  do {
		iter++;
		status = gsl_root_fsolver_iterate(s);
		r=gsl_root_fsolver_root(s);
		low = gsl_root_fsolver_x_lower(s);
		high = gsl_root_fsolver_x_upper(s);
		status = gsl_root_test_interval(low,high,0,0.001);
	  } while (status==GSL_CONTINUE && iter < max_iter);
	  return(r);
	}
  else
	{
		if(Vbi_r>1) // too large Vgs
		  return(0);
		else
		  return(1);
	}
}

double r_cMESFET (double r,double V, double Vgs, param_cMESFET p)
{
  double v0;
  double VbiV;
  VbiV=p.Vbi-Vgs+V;
  v0 =4*(p.eps_semi*EPSILON)/(EC*p.Nd*POW2(p.radius))*VbiV;
  return(1-r+r*log(r)-v0);
}

double r_cMESFET_gsl(double r, void *pp)
{
  struct R_cMESFET_param * params = (struct R_cMESFET_param *)pp;
  double V = (params->V);
  double Vgs = (params->Vgs);
  param_cMESFET p = (params->p);
  return(r_cMESFET(r,V,Vgs,p));
}

// effect of source and drain resisitance
double func_Ids_cMESFET_R(double Vds,double Vgs,param_cMESFET p,param_solver ps)
{
  double idsmod;
  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;
  int status;
  /* size_t i; */
  size_t iter=0;
  const size_t n=2;
  gsl_multiroot_function F;
  param_cMOSFET pMOS;
  struct Ids_params params = {Vds, Vgs, pMOS,p,ps};
  double x_init[2]={1,0};
  gsl_vector *x = gsl_vector_alloc(n);
  
  F.f = &Ids_cMESFET_RmodFunc;
  F.n = 2;
  F.params = &params;

  gsl_vector_set(x,0,x_init[0]);
  gsl_vector_set(x,1,x_init[1]);
  T=gsl_multiroot_fsolver_hybrids;
  s=gsl_multiroot_fsolver_alloc(T,2);

  gsl_multiroot_fsolver_set(s,&F,x);

  do {
	iter++;
	status = gsl_multiroot_fsolver_iterate(s);
	if(status)
	  break;
	status = gsl_multiroot_test_residual(s->f,1e-7);
  } while(status==GSL_CONTINUE && iter <1000);
  idsmod = Ids_cMESFET(gsl_vector_get(s->x,0),gsl_vector_get(s->x,1),p,ps);
  gsl_multiroot_fsolver_free(s);
  gsl_vector_free(x);
  return(idsmod);
}
	
int Ids_cMESFET_RmodFunc(gsl_vector *x, void *pp, gsl_vector *f)
{
  struct Ids_params * params = (struct Ids_params *) pp;
  const double Vds = (params->Vds);
  const double Vgs = (params->Vgs);
  const param_cMESFET p = (params->cMES);
  const param_solver ps = (params->ps);
  const double x_Vds = gsl_vector_get(x,0);
  const double x_Vgs = gsl_vector_get(x,1);
  
  gsl_vector_set(f,0,x_Vgs-Vgs+Ids_cMESFET(x_Vds,x_Vds,p,ps)*p.Rs);
  gsl_vector_set(f,1,x_Vds-Vds+Ids_cMESFET(x_Vds,x_Vds,p,ps)*(p.Rs+p.Rd));

  return GSL_SUCCESS;
}

/***************************************************/
// Exporter
double func_Qcharge_cMOSFET(double V, double Vgs, param_cMOSFET p, param_solver ps)
{
  return(qroot0(V, Vgs, p, ps));
}

double func_Qcharge2_cMOSFET(double V, double Vgs, param_cMOSFET p)
{
  return(Q_approx(V,Vgs,p));
}

double func_rootfind_Q_cMOSFET(double qq,double V, double Vgs, param_cMOSFET p)
{
  return(qfunc_cMOSFET(qq,V, Vgs, p));
}

double func_Ids_cMOSFET(double Vds,double Vgs,param_cMOSFET p, param_solver ps)
{
  return(Ids_cMOSFET0(Vds,Vgs,p, ps));
}

double func_Ids2_cMOSFET(double Vds,double Vgs,param_cMOSFET p)
{
  return(Ids0_cMOSFET(Vds,Vgs,p));
}

double func_Ids_cMOSFET_R(double Vds,double Vgs,param_cMOSFET p, param_solver ps)
{
  return(Ids_cMOSFET_R0(Vds,Vgs, p, ps));
}

double func_Ids2_cMOSFET_R(double Vds,double Vgs,param_cMOSFET p)
{
  return(Ids0_cMOSFET_R(Vds,Vgs,p));
}

