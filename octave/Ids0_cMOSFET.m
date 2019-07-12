function y=Ids0_cMOSFET(Vds,Vgs,Vt0,dVt0,Vth,Q0,Cox))
  QS = Q_approx(0,Vgs,Vt0,dVt0,Vth,Q0,Cox)
  QD = Q_approx(Vds,Vgs,Vt0,dVt0,Vth,Q0,Cox)
  y=Ids00_cMOSFET(QS,QD)
endfunction
