%module ballistic

%{
#define SWIG_FILE_WITH_INIT
  
#include "../density1d.h"
#include "../ballistic.h"
#include "params.h"
  
  extern double density1d_parabollic00(double EFermi, double Enm, double ems,double temp);
  extern double density1d_nonpara00(double EFermi, double Enm, double alpha_nm, double ems_nm, double temp);
  
  extern double alphaNP00(double Eg, double ems);
  extern double gamma_nm00(double alpha,double ems,double W1,double W2,int n,int m);
  extern double E_nm0(double alpha, double gamma_nm);
  extern double alpha_nm0(double alpha, double gamma_nm);
  extern double ems_nm0(double ems, double gamma_nm);
  
  extern double density1d_all00(double EFermi,double alpha, double ems, double temp,
								double W1, double W2, int nmax, int mmax);
  
  %}

/* %include <numpy.i> */
/* %init %{ */
/*   import_array(); */

%include "../density1d.h"
%include "../ballistic.h"
%include "params.h"

extern double density1d_parabollic00(double EFermi, double Enm, double ems,double temp);
extern double density1d_nonpara00(double EFermi, double Enm, double alpha_nm, double ems_nm, double temp);

extern double alphaNP00(double Eg, double ems);
extern double gamma_nm00(double alpha,double ems,double W1,double W2,int n,int m);
extern double E_nm0(double alpha, double gamma_nm);
extern double alpha_nm0(double alpha, double gamma_nm);
extern double ems_nm0(double ems, double gamma_nm);

extern double density1d_all00(double EFermi,double alpha, double ems, double temp,
							  double W1, double W2, int nmax, int mmax);
/* %inline %{ */

/* %} */

