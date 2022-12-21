function [Fnl, ft] = NLEVAL(m, t, U, Udot,tol, varargin)
%NLEVAL evaluates the nonlinearities in the time domain for given set of
%points

    Nt = length(t);
    if size(U,1)==Nt
        U = U';
    end  % (Nd, Nt)
    if size(Udot,1)==Nt
        Udot = Udot';
    end
    
    if ~isempty(varargin)
        ITMAX = varargin{1};
    else
        ITMAX = 100;
    end
    
    Fnl = zeros(m.Ndofs, Nt);
    for ni=1:length(m.NLTs)
        unlt = (m.NLTs(ni).L*U)';
        unldot = (m.NLTs(ni).L*U)';
        
        Ndnl = size(m.NLTs(ni).L, 1);
        
        if mod(m.NLTs(ni).type,2)==0  % Inst. force
            try
                ft = m.NLTs(ni).func(t(:), unlt, unldot);
            catch me  % not vectorized maybe
                ft = zeros(size(unlt));
                for i=1:size(unlt,1)
                    ft(i,:) = m.NLTs(ni).func(t(i), unlt(i,:)', unldot(i,:));
                end
            end
            % (Nt,Ndnl); (Nt,Ndnl); (Nt,Ndnl) (point-wise)
            
        else  % Hysteretic force
            ft = zeros(Nt, Ndnl);
%             dfdu = zeros(Nt, Ndnl, Ndnl, 1);
% 
%             its  = 0;
%             while its==0 || abs(fprev-ft(end, :))>tol
%                 fprev = ft(end, :);
%                 for ti=1:Nt
%                     tm1 = mod(ti-2, Nt)+1;
%                     [ft(ti,:), dfdu(ti,:,:,:)] = ...
%                         m.NLTs(ni).func(t(ti), unlt(ti,:), 0, t(tm1), ...
%                         unlt(tm1,:), ft(tm1,:), dfdu(tm1,:,:,:));
%                 end
%                 its = its+1;
%             end
            dfdu = zeros(1, Ndnl, Ndnl, 1);
            
            its  = 0;
            while its==0 || any(abs(fprev-ft(end, :))>tol)
                fprev = ft(end, :);
                for ti=1:Nt
                    tm1 = mod(ti-2, Nt)+1;
                    ft(ti,:) = m.NLTs(ni).func(t(ti), unlt(ti,:)', ...
                        t(tm1), unlt(tm1,:)', ft(tm1,:), 0, 0);
%                     fprintf('%d, ', ti);
                end
%                 fprintf('\n');
                its = its+1;
                if its>=ITMAX
                    break;
                end
            end
        end
        
        if m.NLTs(ni).type<=5
            Fnl = Fnl + m.NLTs(ni).L'*ft';
        else
            Fnl = Fnl + m.NLTs(ni).Lf*ft';
        end
    end    
    
    Fnl = Fnl';  % Ndofs x Nt
end