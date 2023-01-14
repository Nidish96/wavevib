function [Ns] = QPTIMEINTERP(taus, h)
%QPTIMEINTERP
%
%   USAGE :
%       taus    : (Np, Nc)
%       h       : (Nh, Nc)
%   OUTPUTS :
%       Ns      : (Np, Nhc)
    if sum(all(h==0, 2)) && ~all(h(1,:)==0)
        error('If you want the dc term, put it in the beginning of h')
    end

    Nc = size(h,2);
    Nh = size(h,1);
    Nhc = sum(all(h==0, 2)+2*any(h~=0, 2));
    Np = size(taus, 1);
    
    Ns = zeros(Np, Nhc);
    for ip=1:Np
        k = 1;
        if all(h(1,:)==0)
            Ns(ip, k) = 1;
            k = k+1;
        end
        Ns(ip, k:2:end)   = cos(h(k:end,:)*taus(ip,:)'); 
        Ns(ip, k+1:2:end) = sin(h(k:end,:)*taus(ip,:)'); 
    end
end