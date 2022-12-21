function [R, dRdU, dRdw] = SDOF_NLHYST_HBRESFUN(Uw, m, c, k, Fl, fnl, h, Nt)
%SDOF_NLHYST_HBRESFUN returns the HB residue for an SDOF system
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
  
  % Time Marching for Hysteretic Nonlinearity
  ft = zeros(size(t));
  dfdxt = zeros(size(t));
  dfdxdt = zeros(size(t));
  z = zeros(size(t));
  dzdxt = zeros(size(t));

  ftprev = ft;

  it = 0;
  
  dzdai = zeros(length(t), Nhc);
  jt = zeros(length(t), Nhc);
  jt2 = zeros(length(t), Nhc);  
  
  % while mean((ftprev-ft).^2) > 1e-8 || it == 0
  while it < 2
    ftprev = ft;
    for i = 1:Nt
      [ft(i), z(i), dfdxt(i), dfdxdt(i), dzdxt(i)] = fnl(t(i), ut(i), z(mod(i-2,Nt)+1), udt(i));
      if dzdxt(i)==0
	dzdai(i, :) = dzdai(mod(i-2,Nt)+1, :);
      else
	dzdai(i, :) = dzdxt(i)*cst(i, :);
      end
      if z(i)==z(mod(i-2,Nt)+1)  % stick
	jt(i, :) = dfdxt(i)*(cst(i, :)-cst(mod(i-2,Nt)+1,:)) + jt(mod(i-2,Nt)+1,:);
      else  % slip
	jt(i, :) = 0;
      end
    end
    it = it+1;
  end
  
  % End of time marching

  Fnl = GETFOURIERCOEFF(h, ft);
		  %  Jnl = GETFOURIERCOEFF(h, dfdxt.*cst+dfdxdt.*sct);
  jt = dfdxt.*(cst-dzdai);
  Jnl = GETFOURIERCOEFF(h, jt);
keyboard
  % Residue
  R = E*Uw(1:end-1) + Fnl - Fl;
  dRdU = E + Jnl;
  dRdw = dEdw*Uw(1:end-1);
end
