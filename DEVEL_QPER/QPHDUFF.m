function [FNL, dFNLdU, dFNLdw] = QPHDUFF(Uw, kJs, cJs, gJs, h, Nt)
%HDUFF returns the Fourier Coefficients of the cubic spring along with
%linear connections for given displacement coefficients.
%
%   USAGE: 
%       [FNL, dFNLdU, dFNLdw] = HDUFF(Uw, kJ, cJ, gJ, h, Nt);
%   INPUTS:
%       Uw      : (Nhc+Nc,1)
%       kJ,cJ,gJ: (scalar)
%       h       : (Nh,1)
%       Nt      : (int)
%   OUTPUTS:
%       FNL     : (Nhc,1)
%       dFNLdU  : (Nhc,Nhc)
%       dFNLdw  : (Nhc,1)
    
    Nhc = sum(all(h==0, 2)+2*any(h~=0, 2));
    Nc = length(Uw)-Nhc;
    ws = Uw(end-Nc+1:end);

    D1 = kron(diag(h*ws), [0 1;-1 0]);  % Fourier Differentiation mx
    dD1w = reshape(cell2mat(arrayfun(@(a) kron(diag(h(:,a)), [0 1;-1 0]), 1:Nc, 'UniformOutput', false)), Nhc, Nhc, Nc);
    if all(h(1,:)==0)
        D1 = D1(2:end,2:end);
        dD1w = dD1w(2:end, 2:end, :);
    end
    Nd = length(Uw(1:end-Nc))/Nhc;
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

    cst = QPAFT(eye(Nhc), h, Nt, 'f2t');
    sct = QPAFT(D1, h, Nt, 'f2t');

    ut = QPAFT(reshape(Uw(1:end-Nc), Nd, Nhc)', h, Nt, 'f2t');
    udt = QPAFT(D1*reshape(Uw(1:end-Nc), Nd, Nhc)', h, Nt, 'f2t');
    dudtw = cell2mat(arrayfun(@(a) QPAFT(dD1w(:,:,a)*reshape(Uw(1:end-Nc), Nd, Nhc)', h, Nt, 'f2t'), 1:Nc, 'UniformOutput', false));
    
    ft = ut*kJs' + udt*cJs' + ut.^3*gJs';
    dfdu = kron(ones(Nt^Nc,1), kJs');
    for di=1:Nd
        dfdu(di:Nd:end,:) = dfdu(di:Nd:end,:)+3*ut(:,di).^2.*gJs(:,di)';
    end
    dfdud = kron(ones(Nt^Nc,1), cJs');
    
    FNL = reshape(QPAFT(ft, h, Nt, 't2f')', Nd*Nhc,1);
    dFNLdU = zeros(Nhc*Nd);
    dFNLdw = zeros(Nhc*Nd,Nc);
    for di=1:Nd
        for dj=1:Nd
            dFNLdU(di:Nd:end, dj:Nd:end) = ...
                QPAFT(dfdu(dj:Nd:end, di).*cst + dfdud(dj:Nd:end, di).*sct, h, Nt, 't2f');
        end
        dFNLdw(di:Nd:end,:) = QPAFT(permute(sum(reshape(dfdud(:, di), Nd, Nt^Nc)'.*permute(dudtw, [1 3 2]),2), [1 3 2]), h, Nt, 't2f');
    end

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