%module(docstring="I-V characteritics of various FETs") fetmodel

%{
#define SWIG_FILE_WITH_INIT
  
#include "../src/ccm.h"
#include "../src/density1d.h"
#include "../src/density2d.h"
#include "../src/ballistic.h"
#include "../src/capacitor.h"
#include "../src/mos1d.h"
#include "fetmodel.h"

  extern void Qapprox_cMOS_func(double *in_array,double *out_array, int size, param_cMOSFET p);
  extern void Q_cMOS_func(double *in_array,double *out_array, int size, param_cMOSFET p);

  extern void Ids0_cMOS_func(double *in_array,double *out_array, int size, double Vds,param_cMOSFET p);
  extern void Ids_cMOS_func(double *in_array,double *out_array, int size, double Vds,param_cMOSFET p);
  extern void Ids0_cMOS_R_func(double *in_array,double *out_array, int size, double Vds,param_cMOSFET p);
  extern void Ids_cMOS_R_func(double *in_array,double *out_array, int size, double Vds,param_cMOSFET p);
  
  extern void Ids_cMES_func(double *in_array,double *out_array, int size, double Vds,param_cMESFET p);
  extern void Ids_cMES_R_func(double *in_array,double *out_array, int size, double Vds,param_cMESFET p);

  %}

%include <numpy.i>
%init %{
  import_array();
%}

%feature("autodoc","1");

%include "../src/ccm.h"
%include "../src/density1d.h"
%include "../src/density2d.h"
%include "../src/ballistic.h"
%include "../src/capacitor.h"
%include "../src/mos1d.h"
%include "fetmodel.h"

%apply (double* IN_ARRAY1, int DIM1) {(double * in_array, int size_in)}
/* %apply (double* ARGOUT_ARRAY1, int DIM1) {(double * out_array, int size_out)} */

%apply (double* INPLACE_ARRAY1, int DIM1) {(double * out_array, int size_out)}
/* %apply (double** ARGOUTVIEWM_ARRAY2, int* DIM1, int* DIM2){(double** outmat, int* mx, int* my)} */

%inline %{
  void Qapprox_cMOS(double * in_array, int size_in, double * out_array, int size_out,
					param_cMOSFET p) {
	Qapprox_cMOS_func(in_array,out_array,size_in,p);
  }
  void Q_cMOS(double * in_array, int size_in, double * out_array, int size_out,
			  param_cMOSFET p)
  {
	Q_cMOS_func(in_array, out_array, size_in, p);
  }
  void Ids_cMOS(double * in_array, int size_in, double * out_array, int size_out,
				double Vds, param_cMOSFET p) {
	Ids_cMOS_func(in_array, out_array, size_in, Vds, p);
  }
  void Ids0_cMOS(double * in_array, int size_in, double * out_array, int size_out,
				 double Vds, param_cMOSFET p) {
	Ids0_cMOS_func(in_array, out_array, size_in, Vds, p);
  }
  void Ids_cMOS_R(double * in_array, int size_in, double * out_array, int size_out,
				  double Vds, param_cMOSFET p) {
	Ids_cMOS_R_func(in_array, out_array, size_in, Vds, p);
  }
  void Ids0_cMOS_R(double * in_array, int size_in, double * out_array, int size_out,
				   double Vds, param_cMOSFET p) {
	Ids0_cMOS_R_func(in_array, out_array, size_in, Vds, p);
  }
  void Ids_cMES(double * in_array, int size_in, double * out_array, int size_out,
				double Vds, param_cMESFET p) {
	Ids_cMES_func(in_array, out_array, size_in, Vds, p);
  }
  void Ids_cMES_R(double * in_array, int size_in, double * out_array, int size_out,
				  double Vds, param_cMESFET p) {
	Ids_cMES_R_func(in_array, out_array, size_in, Vds, p);
  }
  %}

%pythoncode %{
  from .define_params import (
							  parameters_ballistic,
							  parameters_cMOSFET,
							  )
  from .wrapper import (
  					  density1d_rect1d_all,
  					  density1d_rect1dNP_all,
  					  func_for_findroot_E0_rect1dNP,
  					  E0_rect1dNP_root,
  					  Ids_ballistic1d_rect1dNP,
					  func_for_findroot_E0_2d,
					  E0_2d_root,
					  Ids_ballistic2d_E0,
					  Ids_ballistic2d,
					  func_Qcharge1dMOS0,
  					  )
%}
  
