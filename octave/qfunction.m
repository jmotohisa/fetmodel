function y=qfunction(V,Vgs,V0,Vth,Q0,Cox)
% approximate form of Q (eq.17 in Iniguez Trans. ED)
  qp0 = q_approx0(V,V0,0)
  Vt0=v0 + 2*Vth*ln(1+qp0/q0)
  dVt0 = (2*Cox*Vth*Vth/Q0*qp0/(q0+qp0))
  y=q_approx0(V,Vgs,Vt0,dVt0,Vth,Q0,Cox)
endfunction

#// unified approximated formula for Q (eq. 13 in Iniguez et al., Trans. ED)
