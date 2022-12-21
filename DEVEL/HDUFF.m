function [FNL, dFNLdU, dFNLdw] = HDUFF(Uw, kJ, cJ, gJ, h, Nt)
%HDUFF returns the Fourier Coefficients of the cubic spring along with
%linear connections for given displacement coefficients.
%
%   USAGE: 
%       [FNL, dFNLdU, dFNLdw] = HDUFF(Uw, kJ, cJ, gJ, h, Nt);
%   INPUTS:
%       Uw      : (Nhc+1,1)
%       kJ,cJ,gJ: (scalar)
%       h       : (Nh,1)
%       Nt      : (int)
%   OUTPUTS:
%       FNL     : (Nhc,1)
%       dFNLdU  : (Nhc,Nhc)
%       dFNLdw  : (Nhc,1)
    
    Nhc = sum((h==0)+2*(h~=0));
    D1 = kron(diag(h),[0 1;-1 0])*Uw(end);  % Fourier Differentiation mx
    if h(1)==0
        D1 = D1(2:end,2:end);
    end

    ut = AFT(Uw(1:end-1), Nt, h, 'f2t');
    udt = AFT(D1*Uw(1:end-1), Nt, h, 'f2t');
    
    cst = AFT(eye(Nhc), Nt, h, 'f2t');
    sct = AFT(D1, Nt, h, 'f2t');

    ft = kJ*ut + cJ*udt + gJ*ut.^3;
    dfdu = kJ + 3*gJ*ut.^2;
    dfdud = cJ*ones(Nt,1);

    FNL = AFT(ft, Nt, h, 't2f');
    dFNLdU = AFT(dfdu.*cst + dfdud.*sct, Nt, h, 't2f');
    dFNLdw = AFT(dfdud.*udt/Uw(end), Nt, h, 't2f');
end