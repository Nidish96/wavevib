function [R, dRdU, dRdw] = SDOF_NL_HBRESFUN(Uw, m, c, k, Fl, fnl, h, Nt)
%SDOF_NL_HBRESFUN returns the HB residue for an SDOF system
% 
% 
%

  Nhc = sum(h==0)+2*sum(h~=0);
  
  [E, dEdw] = HARMONICSTIFFNESS(m, c, k, Uw(end), h);
  D1 = HARMONICSTIFFNESS(0, 1, 0, Uw(end), h);

  t = linspace(0, 2*pi, Nt+1)'; t(end) = [];
  ut = TIMESERIES_DERIV(Nt, h, Uw(1:end-1), 0);
  udt = TIMESERIES_DERIV(Nt, h, Uw(1:end-1), 1);

  cst = TIMESERIES_DERIV(Nt, h, eye(Nhc), 0);
  sct = TIMESERIES_DERIV(Nt, h, D1, 0);

  [ft, dfdxt, dfdxdt] = arrayfun(fnl, t, ut, udt);

  Fnl = GETFOURIERCOEFF(h, ft);
  Jnl = GETFOURIERCOEFF(h, dfdxt.*cst+dfdxdt.*sct);

  % Residue
  R = E*Uw(1:end-1) + Fnl - Fl;
  dRdU = E + Jnl;
  dRdw = dEdw*Uw(1:end-1);
end
