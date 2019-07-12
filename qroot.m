function y=qroot0(V,Vgs)
  low=0
  high=1-low
								#	FindRoots/Q/L=(low) qfunc_cMOSFET,param_cMOSFET
  [x,fval,info]=fsolve(@(x) qfunc_cMOSFET(x,Vgs,Vth,dphi,delta,radius,Cox,Q0),[1.0;0]);
endfunction
