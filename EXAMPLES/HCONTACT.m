function [FNL, dFNLdU, dFNLdw] = HCONTACT(Uw, knl, gap, h, Nt)
%HCONTACT returns the Fourier Coefficients of the contact model implemented
%as a linear-penalty model. Force that is returned is:
%                   knl*max(u-gap, 0)
%   positive (u-gap) quantity is considered as contact and negative is
%   considered as separation (zero force). 
%
%   USAGE: 
%       [FNL, dFNLdU, dFNLdw] = HCONTACT(Uw, knl, gap, h, Nt);
%   INPUTS:
%       Uw      : (Nhc*Nd+1,1) Vector of nonlinear relative DOFs+freq.
%       knl,gap : (scalar) or (Nd,1) parameters
%       h       : (Nh,1) Vector of Harmonics
%       Nt      : (int) Number of samples for AFT
%   OUTPUTS:
%       FNL     : (Nhc*Nd,1) Force harmonics
%       dFNLdU  : (Nhc*Nd,Nhc*Nd) Force harmonic jacobian wrt U
%       dFNLdw  : (Nhc*Nd,1) Force harmonic jacobian wrt w
    
    Nhc = sum((h==0)+2*(h~=0));
    cst = AFT(eye(Nhc), h, Nt, 'f2t');

    Nd = (size(Uw,1)-1)/Nhc;  % Number of harmonics
    if length(knl)==1
        knl = ones(1,Nd)*knl;
    end
    if length(gap)==1
        gap = ones(1,Nd)*gap;
    end

    % Evaluate nonlinear force through AFT
    ut = AFT(reshape(Uw(1:end-1), Nd, Nhc)', h, Nt, 'f2t');
    ft = knl.*max(ut-gap, 0);
    dfdu = knl.*(ft~=0);

    FNL    = reshape(AFT(ft, h, Nt, 't2f')', Nd*Nhc,1);
    dFNLdU = zeros(Nd*Nhc);
    for di=1:Nd
        dFNLdU(di:Nd:end,di:Nd:end) = AFT(dfdu(:,di).*cst, h, Nt, 't2f');
    end
    dFNLdw = zeros(Nhc*Nd,1);
end