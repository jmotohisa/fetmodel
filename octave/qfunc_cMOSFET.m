function y=qfunc_cMOSFET(qq,Vgs,Vth,dphi,delta,radius,Cox,Q0)
  qqq1 = Vgs-dphi-V-Vth*ln(8/(delta*radius^2))
  qqq2=- (qq/Cox + Vth*(ln(qq/Q0)+ln(1+(qq/Q0))))
  y=qqq1+qqq2
endfunction
