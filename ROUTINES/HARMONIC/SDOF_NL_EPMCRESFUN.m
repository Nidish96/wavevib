function [R, dRdUwx, dRda] = SDOF_NL_EPMCRESFUN(Uwxa, m, c, k, Fl, fnl, h, Nt)
%SDOF_NL_EPMCRESFUN returns the EPMC residue for an SDOF system
% 
% 
%

  Nhc = sum(h==0)+2*sum(h~=0);
  
  lA = Uwxa(end);
  A = 10^lA;  % EPMC amplitude
  dAdlA = A*log(10);
  
  xi = Uwxa(end-1);
  
  [E, dEdw] = HARMONICSTIFFNESS(m, c-xi*m, k, Uwxa(end-2), h);
  dEdxi = HARMONICSTIFFNESS(0, -m, k, Uwxa(end-2), h);
  
  D1 = HARMONICSTIFFNESS(0, 1, 0, Uwxa(end-2), h);

  t = linspace(0, 2*pi, Nt+1)'; t(end) = [];
  ut = TIMESERIES_DERIV(Nt, h, A*Uwxa(1:end-3), 0);
  udt = TIMESERIES_DERIV(Nt, h, A*Uwxa(1:end-3), 1);

  cst = TIMESERIES_DERIV(Nt, h, eye(Nhc), 0);
  sct = TIMESERIES_DERIV(Nt, h, D1, 0);

  [ft, dfdxt, dfdxdt] = arrayfun(fnl, t, ut, udt);

  Fnl = GETFOURIERCOEFF(h, ft);
  Jnl = GETFOURIERCOEFF(h, dfdxt.*cst+dfdxdt.*sct);

  % Residue
  R = [E*(A*Uwxa(1:end-3)) + Fnl;
       m*(Uwxa(2)^2+Uwxa(3)^2)-1;
       Fl'*Uwxa(1:end-3)];
   
  dRdUwx = [E*A+Jnl*A, dEdw*(A*Uwxa(1:end-3)), dEdxi*(A*Uwxa(1:end-3));
          2*m*[0 Uwxa(2) Uwxa(3) zeros(size(Uwxa(4:end-1)'))];
          [Fl' 0 0]];
  dRda = [(E+Jnl)*Uwxa(1:end-3)*dAdlA; 0; 0];
end
