function [R, dRdUwx, dRda] = EPMCRESFUN(m, Uwxa, Fl, h, Nt, tol, varargin)
%EPMCRESFUN
%
%   USAGE: 
%       [R, dRdUwx, dRda] = EPMCRESFUN(m, Uwxa, Fl, h, Nt, tol, varargin)
%   INPUTS:
%       Uwxa    :
%       Fl      :
%       h       :
%       Nt      :
%       tol     :
%   OUTPUTS:    
%       R       :
%       dRdUwx  :
%       dRda    :

    Nhc = sum((h==0)+2*(h~=0));

    la = Uwxa(end);
    A = 10^la;
    dAdla = A*log(10);
    
    h0 = double(h(1)==0);
    Asc = kron([ones(h0,1); A*ones(Nhc-h0,1)], ones(m.Ndofs,1));
    dAscdA = kron([zeros(h0,1); ones(Nhc-h0,1)], ones(m.Ndofs,1));
    
    xi = Uwxa(end-1);
    w = Uwxa(end-2);
    
    [E, dEdw] = HARMONICSTIFFNESS(m.M, m.C-xi*m.M, m.K, w, h);
    dEdxi = HARMONICSTIFFNESS(m.M*0, -m.M, m.M*0, w, h);
    
    t = linspace(0, 2*pi, Nt+1)';  t(end) = [];
    cst = TIMESERIES_DERIV(Nt, h, eye(Nhc), 0);  
    sct = w*TIMESERIES_DERIV(Nt, h, eye(Nhc), 1);
    
    FNL = zeros(m.Ndofs*Nhc, 1);
    dFNL = zeros(m.Ndofs*Nhc);
    for ni=1:length(m.NLTs)
        Unl = (m.NLTs(ni).L*reshape(Asc.*Uwxa(1:end-3), m.Ndofs, Nhc))';  % Nhc x Ndnl
        Ndnl = size(m.NLTs(ni).L, 1);
        
        unlt = TIMESERIES_DERIV(Nt, h, Unl, 0);  % Nt x Ndnl
        unldot = w*TIMESERIES_DERIV(Nt, h, Unl, 1);  % Nt x Ndnl
        
        if mod(m.NLTs(ni).type, 2)==0  % Instantaneous force
            [ft, dfdu, dfdud] = m.NLTs(ni).func(t, unlt, unldot);
            % (Nt,Ndnl); (Nt,Ndnl); (Nt,Ndnl) (point-wise)
            F = GETFOURIERCOEFF(h, ft);
            J = zeros(size(m.NLTs(ni).L,1)*Nhc, size(m.NLTs(ni).L,1)*Nhc);
            dFdU = reshape(GETFOURIERCOEFF(h, reshape(dfdu.*permute(cst, [1, 3, 2]), Nt, Ndnl*Nhc)),...
                Nhc, Ndnl, Nhc);
            for di=1:Ndnl
                J(di:Ndnl:end, di:Ndnl:end) = dFdU(:, di, :);
            end
        else  % Hysteretic force
            ft = zeros(Nt, Ndnl);
            dfdu = zeros(Nt, Ndnl, Ndnl, Nhc);
            
            its = 0;
            
            while its==0 || max(abs(fprev-ft(end, :)))>tol
                fprev = ft(end, :);
                for ti=1:Nt
                    tm1 = mod(ti-2, Nt)+1;
                    [ft(ti,:), dfdu(ti,:,:,:)] = ...
                        m.NLTs(ni).func(t(ti), unlt(ti,:), h, t(tm1), ...
                        unlt(tm1,:), ft(tm1,:), dfdu(tm1,:,:,:));
                end
                its = its+1;
%                keyboard 
            end
%             fprintf('%d\n', its);
            F = GETFOURIERCOEFF(h, ft);
            J = zeros(size(m.NLTs(ni).L,1)*Nhc, size(m.NLTs(ni).L,1)*Nhc);
            for di=1:Ndnl
                for dj=1:Ndnl
                    tmp = squeeze(dfdu(:, di, dj, :));
                    if ~isempty(find(tmp~=0, 1))
                        J(di:Ndnl:end, dj:Ndnl:end) = ...
                            GETFOURIERCOEFF(h, tmp);
                    end
                end
            end
        end
        
        if m.NLTs(ni).type<=5  % Self adjoint forcing
            FNL = FNL + reshape(m.NLTs(ni).L'*F', Nhc*m.Ndofs, 1);
            dFNL = dFNL + kron(eye(Nhc), m.NLTs(ni).L')*J*kron(eye(Nhc), m.NLTs(ni).L);
	    else  % Non-self adjoint forcing
            FNL = FNL + reshape(m.NLTs(ni).Lf*F', Nhc*m.Ndofs, 1);
            dFNL = dFNL + kron(eye(Nhc), m.NLTs(ni).Lf)*J*kron(eye(Nhc), m.NLTs(ni).L);                
        end
    end
    
    % Residue
    h0 = double(h(1)==0);
    Fstat = kron([ones(h0,1); zeros(Nhc-h0,1)], ones(m.Ndofs,1)).*Fl;  % Only static loads if it's there
    Fdyn = kron([zeros(h0,1); ones(Nhc-h0, 1)], ones(m.Ndofs,1)).*Fl;  % Only dynamic loads (for phase constraint)
    
%   Use all harmonic amplitudes as mode shape and normalize by sum of all harmonics
%     R = [E*(Asc.*Uwxa(1:end-3))+FNL-Fstat;
%         Uwxa(1:end-3)'*kron(diag([0 ones(1,Nhc-1)]), m.M)*Uwxa(1:end-3)-1;
%         Fdyn'*Uwxa(1:end-3)];
%     dRdUwx = [(E+dFNL)*diag(Asc), dEdw*(Asc.*Uwxa(1:end-3)), dEdxi*(Asc.*Uwxa(1:end-3));
%         2*Uwxa(1:end-3)'*kron(diag([0 ones(1,Nhc-1)]), m.M), 0, 0;
%         Fdyn', 0, 0];
%     dRda = [(E+dFNL)*(dAscdA.*Uwxa(1:end-3))*dAdla;
%         0; 
%         0];

%   Use first harmonic amplitudes as mode shape
    R = [E*(Asc.*Uwxa(1:end-3))+FNL - Fstat;
        Uwxa(m.Ndofs*h0+(1:2*m.Ndofs))'*blkdiag(m.M,m.M)*Uwxa(m.Ndofs*h0+(1:2*m.Ndofs))-1;
        Fdyn'*Uwxa(1:end-3)];
    dRdUwx = [(E+dFNL)*diag(Asc), dEdw*(Asc.*Uwxa(1:end-3)), dEdxi*(Asc.*Uwxa(1:end-3));
        2*Uwxa(1:end-3)'*kron(diag([0 ones(1,Nhc-1)]), m.M), 0, 0;
        Fdyn', 0, 0];
    dRda = [(E+dFNL)*(dAscdA.*Uwxa(1:end-3))*dAdla;
        0; 
        0];
    
    if ~isempty(varargin)
        dRdUwx = [dRdUwx dRda];
    end
end