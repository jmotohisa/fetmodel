function y=Q_approx(V,Vgs,Vth,V0,Q0,Cox)
% approximate form of Q (eq.17 in Iniguez Trans. ED)
  qp0 = q_approx0(V,Vgs,V0,0)
  Vt0=V0 + 2*Vth*ln(1+qp0/Q0)
  dVt0 = (2*Cox*Vth^2/Q0)*qp0/(Q0+qp0)
  y=q_approx0(V,Vgs,Vt0,dVt0,Vth,Q0,Cox)
endfunction
