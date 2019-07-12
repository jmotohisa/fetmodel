%module pycfet

%{
#define SWIG_FILE_WITH_INIT
  
  #incldde "../ccm.h"
  #include "pycfet.h"
%}

  /* %include <ctl.h> */
  %include "../ccm.h"
  %include "pycfet.h"

extern double Ids_cMOS(double Vds, double Vgs, param_cMOSFET p);
extern param_cMOSFET *param_cMOSFET_new(void );  
extern param_cMESFET *param_cMESFET_new(void );  
