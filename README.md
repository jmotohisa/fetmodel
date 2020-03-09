* Calculation of Drain Current in various FETs

** cylindrical MOSFET: 
Iniguez, B., Jimenez, D., Roig, J., Hamid, H., Marsal, L., & Pallares, J. (2005). Explicit continuous model for long-channel undoped surrounding gate MOSFETs. IEEE Transactions on Electron Devices, 52(8), 1868–1873.

** cylindrical MESFET: see doc

** Ballistic MOSFET


**Basic Function to export (ccm.h)
~~~
double func_Qcharge_cMOSFET(double V, double Vgs, param_cMOSFET p, param_solver ps);
double func_Ids_cMOSFET(double Vds, double Vgs, param_cMOSFET p, param_solver ps);
double func_Ids_cMOSFET_R(double Vds, double Vgs, param_cMOSFET p, param_solver ps);
double func_Qcharge2_cMOSFET(double V, double Vgs, param_cMOSFET p);
double func_Ids2_cMOSFET(double Vds, double Vgs, param_cMOSFET p);
double func_Ids2_cMOSFET_R(double Vds, double Vgs, param_cMOSFET p);
double func_rootfind_Q_cMOSFET(double qq,double V, double Vgs, param_cMOSFET p);
double func_rootfind_logQ_cMOSFET(double qq,double V, double Vgs, param_cMOSFET p);

double Ids_cMESFET(double Vds, double Vgs, param_cMESFET p);
double Ids_cMESFET_R(double Vds, double Vgs, param_cMESFET p);
~~~

*Exported functions in Libctl/Guile (fetmodel.scm/cfet_libctl.c):
~~~~
(Qcharge-cMOSFET Vgs)
(Ids-cMOSFET Vds Vgs)
(Ids-cMOSFET-R Vds Vgs)
(Qcharge2-cMOSFET Vgs)
(Ids2-cMOSFET Vds Vgs)
(Ids2-cMOSFET-R Vds Vgs)
(frf_log_cMOSFET qq V Vgs FET_params)
(frf_logQ_cMOSFET qq V Vgs FET_params)
(func-Ids-cMESFET Vds Vgs)
(func-Ids-cMESFET-R Vds Vgs)
~~~~

**Exported functions in Python ():
~~~
Q_cMOSFET( Vgs, param_cMOSFET p)
Qapprox_cMOSFET( Vgs, param_cMOSFET p)
Ids_cMOSFET( Vds,  Vgs, param_cMOSFET p)
Ids_cMOSFET_R( Vds,  Vgs, param_cMOSFET p)
Ids0_cMOSFET( Vds,  Vgs, param_cMOSFET p)
Ids0_cMOSFET_R( Vds,  Vgs, param_cMOSFET p)
Ids_cMESFET(Vds, Vgs, param_cMESFET p)
Ids_cMESFET_R(Vds, Vgs, param_cMESFET p)

Q_cMOS(numpyarray Vgs, numpyarray Q, param_cMOSFET p)
Qapprox_cMOS(numpyarray Vgs, numpyarray Q, param_cMOSFET p);
Ids_cMOS(numpyarray Vgs, numpyarray Ids, Vds,  param_cMOSFET p);
Ids_cMOS_R(numpyarray Vgs, numpyarray Ids, Vds,param_cMOSFET p);
Ids0_cMOS(numpyarray Vgs, numpyarray Ids, Vds, param_cMOSFET p);
Ids0_cMOS_R(numpyarray Vgs, numpyarray Ids, Vds, param_cMOSFET p);
Ids_cMES(numpyarray Vgs, numpyarray Ids, Vds, param_cMESFET p)
Ids_cMES_R(numpyarray Vgs, numpyarray Ids, Vds, param_cMESFET p)
~~~


