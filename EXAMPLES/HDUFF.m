function [FNL, dFNLdU, dFNLdw] = HDUFF(Uw, kJs, cJs, gJs, h, Nt)
%HDUFF returns the Fourier Coefficients of the cubic spring along with
%linear connections for given displacement coefficients. Force that is
%returned is:
%           kJs * us + cJs * uds + gJs * us.^3
%
%   USAGE: 
%       [FNL, dFNLdU, dFNLdw] = HDUFF(Uw, kJ, cJ, gJ, h, Nt);
%   INPUTS:
%       Uw      : (Nd*Nhc+1,1)
%       kJ,cJ,gJ: (Nd,Nd)
%       h       : (Nh,1)
%       Nt      : (int)
%   OUTPUTS:
%       FNL     : (Nd*Nhc,1)
%       dFNLdU  : (Nd*Nhc,Nd*Nhc)
%       dFNLdw  : (Nd*Nhc,1)
    
    Nhc = sum((h==0)+2*(h~=0));
    D1 = kron(diag(h),[0 1;-1 0])*Uw(end);  % Fourier Differentiation mx
    if h(1)==0
        D1 = D1(2:end,2:end);
    end
    Nd = length(Uw(1:end-1))/Nhc;
    if rem(Nd,1)~=0
        error('Nd not an integer. Value is %f .', Nd);
    end
    if all(size(kJs)~=Nd)
        kJs = kJs(1)*eye(Nd);
    end
    if all(size(cJs)~=Nd)
        cJs = cJs(1)*eye(Nd);
    end
    if all(size(gJs)~=Nd)
        gJs = gJs(1)*eye(Nd);
    end

    cst = AFT(eye(Nhc), h, Nt, 'f2t');
    sct = AFT(D1, h, Nt, 'f2t');

    ut = AFT(reshape(Uw(1:end-1), Nd, Nhc)', h, Nt, 'f2t');
    udt = AFT(D1*reshape(Uw(1:end-1), Nd, Nhc)', h, Nt, 'f2t');
    
    ft = ut*kJs' + udt*cJs' + ut.^3*gJs';
    dfdu = kron(ones(Nt,1), kJs');
    for di=1:Nd
        dfdu(di:Nd:end,:) = dfdu(di:Nd:end,:)+3*ut(:,di).^2.*gJs(:,di)';
    end
    dfdud = kron(ones(Nt,1), cJs');
    
    FNL = reshape(AFT(ft, h, Nt, 't2f')', Nd*Nhc,1);
    dFNLdU = zeros(Nhc*Nd);
    dFNLdw = zeros(Nhc*Nd,1);
    for di=1:Nd
        for dj=1:Nd
            dFNLdU(di:Nd:end, dj:Nd:end) = ...
                AFT(dfdu(dj:Nd:end, di).*cst + dfdud(dj:Nd:end, di).*sct, h, Nt, 't2f');
        end
        dFNLdw(di:Nd:end) = AFT(sum(reshape(dfdud(:, di), Nd, Nt)'.*udt/Uw(end),2), h, Nt, 't2f');
    end

%     % Simple Version
%     ut = AFT(Uw(1:end-1), Nt, h, 'f2t');
%     udt = AFT(D1*Uw(1:end-1), Nt, h, 'f2t');
% 
%     ft = kJs*ut + cJs*udt + gJs*ut.^3;
%     dfdu = kJs + 3*gJs*ut.^2;
%     dfdud = cJs*ones(Nt,1);
% 
%     FNL = AFT(ft, Nt, h, 't2f');
%     dFNLdU = AFT(dfdu.*cst + dfdud.*sct, Nt, h, 't2f');
%     dFNLdw = AFT(dfdud.*udt/Uw(end), Nt, h, 't2f');
end