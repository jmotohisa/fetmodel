radius = 50; # unit in nm
Lg = 0.3; #unit in micron
Cox=0.1 ; #unit in pF/m
mue=3000; # unit in cm^2/sVs
ni=1e15; # unit in cm^-3
temp=300 # unit in K
eps_semi=15.15; 
dphi=0; # unit in eV

eps_ox=20;
tox=20e-9;

g_KBC = 1.38e-23;
g_EC = 1.602e-19;
g_EPSILON =8.85e-12;

radius = radius*1e-9;
Lg = Lg*1e-6;
Cox = Cox*1e-12
mue = mue*1e-4
ni = ni*1e6

KBT=temp*g_KBC/g_EC;
Vth=temp*g_KBC/g_EC;
Cox0 = eps_ox*g_EPSILON/(radius*log(1+tox/radius));
Q0=4*eps_semi*g_EPSILON/radius*Vth;
#	ni = n_intrinsic(temperature,me,mh,Eg,6);
delta=(g_EC*ni/(KBT*(eps_semi*g_EPSILON)));
V0= (dphi + Vth*log(8/(delta*radius^2)));
K = 2*pi * radius/Lg *mue;

Cox=Cox0;

Ids=Ids0_cMOSFET(Vds,Vgs,Vt0,dVt0,Vth,Q0,Cox)*K;
