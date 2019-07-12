%module cfet

%{
#define SWIG_FILE_WITH_INIT
  
  /* #include <ctl.h> */
  #include "../ctl-io.h"
  #include "cfet.h"
%}

  /* %include <ctl.h> */
  %include "../ctl-io.h"
  %include "cfet.h"

extern double Ids_cMOS(double Vds, double Vgs, param_cMOSFET p);
extern param_cMOSFET *param_cMOSFET_new(void );  
extern param_cMESFET *param_cMESFET_new(void );  
