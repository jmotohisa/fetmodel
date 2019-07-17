%module pycfet

%{
#define SWIG_FILE_WITH_INIT
  
  #include "../ccm.h"
  #include "pycfet.h"
%}

%include <numpy.i>
%init %{
  import_array();
%}

%include "../ccm.h"
%include "pycfet.h"

%apply (double* IN_ARRAY1, int DIM1) {(double * in_array, int size_in)}
/* %apply (double* ARGOUT_ARRAY1, int DIM1) {(double * out_array, int size_out)} */

%apply (double* INPLACE_ARRAY1, int DIM1) {(double * out_array, int size_out)}
/* %apply (double** ARGOUTVIEWM_ARRAY2, int* DIM1, int* DIM2){(double** outmat, int* mx, int* my)} */

%inline %{
  void Ids_cMOS(double * in_array, int size_in, double * out_array, int size_out,
				double Vds, param_cMOSFET p) {
	Ids_cMOS_func(in_array, out_array, size_in, Vds, p);
  }
  void Qapprox_cMOS(double * in_array, int size_in, double * out_array, int size_out,
					param_cMOSFET p) {
	Qapprox_cMOS_func(in_array,out_array,size_in,p);
  }
  void Q_cMOS(double * in_array, int size_in, double * out_array, int size_out,
			  param_cMOSFET p)
  {
	Q_cMOS_func(in_array, out_array, size_in, p);
  }
  %}

/* extern param_cMOSFET *param_cMOSFET_new(void );   */
/* extern param_cMESFET *param_cMESFET_new(void );   */

